%% Verify the SALT solution with a true time-domain EMT simulation.
%
% The SALT solve returns a periodic steady state directly. The check is
% that an EMT simulation started from it stays there: integrate period by
% period until consecutive periods line up, then plot the EMT waveforms
% with the SALT waveforms drawn on top. If SALT is right, the two overlap
% and the EMT run reaches periodicity almost immediately.

% Time-domain model-index globals used below (populated by emtInitIndices,
% which the EMT constructor calls)
global idxGenrouTD_Iabc idxGenrouTD_Wr ...
       idxIbrGFMTD_Iabc idxIbrGFMTD_Wr ...
       idxLoadCompTD_Iabc idxLoadCompTD_Wr idxLoadComp_Wr %#ok<GVMIS>

wSys = salt.x_v(salt.getRowIdx_w);
fSys = wSys/(2*pi);

% The EMT verification integrates the full three-phase system period by
% period, which costs about an hour per two periods at 2000 buses. Large
% cases set verifyWithEmtSim false and get the operating-point report only.
if ~exist('verifyWithEmtSim','var')
    verifyWithEmtSim = true;
end
if verifyWithEmtSim

%% Build the EMT twin, initialised from the SALT solution
emt = salt.buildEmt;

deltaT = salt.deltaT;
tSamp  = salt.tSim;
numT   = numel(tSamp);

sysSizeEmt = emt.getSysSize;

%% Simulate until periodic
% Downsampled grid, so long runs stay plottable
[~,tsampIdxPeriodic_downsample,numTPeriodic_downsample,deltaT_downsample] = ...
    getDownsample(numT,deltaT,factor_downsample);

% Rows watched for periodicity: every bus voltage, and every device's
% terminal current
idxLin = 1:(3*emt.numBus);
idxIabc = emt.getRowIdx_fromGenrouIndices(1:emt.numGenrou)-1 + idxGenrouTD_Iabc.';
if emt.numIbrGFM > 0
    idxIabc = [idxIabc(:); reshape( ...
        emt.getRowIdx_fromIbrGFMIndices(1:emt.numIbrGFM)-1 + ...
        idxIbrGFMTD_Iabc.',[],1)];
end
if emt.numLoadComp > 0
    idxIabc = [idxIabc(:); reshape( ...
        emt.getRowIdx_fromLoadCompIndices(1:emt.numLoadComp)-1 + ...
        idxLoadCompTD_Iabc.',[],1)];
end
idxIabc = idxIabc(:);

% Rotor / inverter speeds, which must settle at the SALT system frequency
idxWr = emt.getRowIdx_fromGenrouIndices((1:emt.numGenrou).')-1 + idxGenrouTD_Wr;
wExpected = ones(emt.numGenrou,1)*wSys;
if emt.numIbrGFM > 0
    idxWr = [idxWr; emt.getRowIdx_fromIbrGFMIndices((1:emt.numIbrGFM).')-1 ...
        + idxIbrGFMTD_Wr];
    wExpected = [wExpected; ones(emt.numIbrGFM,1)*wSys];
end
if emt.numLoadComp > 0
    idxWr = [idxWr; emt.getRowIdx_fromLoadCompIndices((1:emt.numLoadComp).')-1 ...
        + idxLoadCompTD_Wr];
    wExpected = [wExpected; salt.x_v( ...
        salt.getRowIdx_fromLoadCompIndices((1:salt.numLoadComp).')-1 ...
        + idxLoadComp_Wr)];
end

% Storage
x_v_lin_prevPeriod  = zeros(numel(idxLin),numT);
x_v_lin_currPeriod  = zeros(numel(idxLin),numT);
x_v_iabc_prevPeriod = zeros(numel(idxIabc),numT);
x_v_iabc_currPeriod = zeros(numel(idxIabc),numT);
x_v_all             = zeros(sysSizeEmt,numT);
x_v_TD_plotting     = zeros(sysSizeEmt, ...
    totIters_periodic*numTPeriodic_downsample);
x_v_TD_plottingWr   = zeros(numel(idxWr), ...
    totIters_periodic*numTPeriodic_downsample);

idx_periodic = 0;
isPeriodic = false;
tic
while (idx_periodic < totIters_periodic) && (~isPeriodic)
    for idxT = 1:numT
        % History currents are fixed across the Newton iterations
        emt.I = zeros(sysSizeEmt,1);
        emt = emt.stampLinear_I(idxT);
        I_linear_static = emt.I;

        idxNR = 0; convergeSim = false;
        while (idxNR<=totIters_NR) && (~convergeSim)
            emt.I = I_linear_static;
            emt = emt.stampNonlinear_G_I;
            x_v_prevGuess_TD = emt.x_v;
            emt.x_v = emt.G \ emt.I;
            convergeSim = isErrSmall(emt.x_v,x_v_prevGuess_TD, ...
                relErrTol_NR,absErrTol_NR);
            idxNR = idxNR + 1;
        end
        if ~convergeSim, warning('NR did not converge!!'), end

        emt.boundThetaRadians;

        % NOTE: For a periodic simulation the "previous" point of the next
        % period must be the LAST point of this period, not the first, so
        % the state is not advanced across the period boundary.
        if idxT~=numT
            emt.x_v_prevTime = emt.x_v;
            emt = emt.setIPrev_fromVPrev;
        end

        x_v_lin_currPeriod(:,idxT)  = emt.x_v(idxLin);
        x_v_iabc_currPeriod(:,idxT) = emt.x_v(idxIabc);
        x_v_all(:,idxT)             = emt.x_v;
    end
    idx_periodic = idx_periodic + 1;

    idxSave = (idx_periodic-1)*numTPeriodic_downsample + ...
        (1:numTPeriodic_downsample);
    x_v_TD_plotting(:,idxSave)   = x_v_all(:,tsampIdxPeriodic_downsample);
    x_v_TD_plottingWr(:,idxSave) = x_v_all(idxWr,tsampIdxPeriodic_downsample);

    % Periodic yet? Voltage, current and speed all have to agree.
    absTolV = absErrMagnitudePercentage_periodic* ...
        max(abs([x_v_lin_currPeriod;x_v_lin_prevPeriod]),[],'all');
    isPeriodicV = isErrSmall(x_v_lin_currPeriod,x_v_lin_prevPeriod, ...
        relErrTol_periodic,absTolV);
    absTolI = absErrMagnitudePercentage_periodic* ...
        max(abs([x_v_iabc_currPeriod;x_v_iabc_prevPeriod]),[],'all');
    isPeriodicI = isErrSmall(x_v_iabc_currPeriod,x_v_iabc_prevPeriod, ...
        relErrTol_periodic,absTolI);
    isPeriodicW = all(abs(x_v_TD_plottingWr(:,idxSave(end))-wExpected)<1e-2);
    isPeriodic = isPeriodicV && isPeriodicI && isPeriodicW;

    if isPeriodic
        % The last point of the period must equal the first point of the next
        assert(isErrSmall(x_v_lin_currPeriod(:,1),x_v_lin_currPeriod(:,end), ...
            relErrTol_periodic,absTolV))
        assert(isErrSmall(x_v_iabc_currPeriod(:,1),x_v_iabc_currPeriod(:,end), ...
            relErrTol_periodic,absTolI))
    else
        x_v_lin_prevPeriod  = x_v_lin_currPeriod;
        x_v_iabc_prevPeriod = x_v_iabc_currPeriod;
    end
end
timeSim = toc;

if ~isPeriodic
    warning(['EMT (SALT-initialised) did not reach steady state after ' ...
        num2str(idx_periodic) ' period(s)...'])
else
    disp(['EMT (SALT-initialised) reached steady state after ' ...
        num2str(idx_periodic) ' period(s).'])
end
disp(['EMT simulation time: ' num2str(timeSim) ' seconds.'])

%% Plot: EMT waveforms with the SALT solution overlaid
numTimepoints = idx_periodic*numTPeriodic_downsample;
tPlot = 0:deltaT_downsample:(deltaT_downsample*(numTimepoints-1));

resultsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
if ~isfolder(resultsDir), mkdir(resultsDir), end

if ~exist('listBusToPlot','var') || isempty(listBusToPlot)
    listBusToPlot = salt.bus;
end

for busNum = reshape(listBusToPlot,1,[])

    % EMT waveform
    idxRow_sim = emt.getRowIdxAbc_fromBusNum(busNum);
    vabc_TD = x_v_TD_plotting(idxRow_sim,1:numTimepoints);

    % SALT waveform, one period reconstructed from the phasor solution and
    % tiled across however many periods the EMT run took.
    % NOTE: The frequencies are fixed beforehand, so the coefficients are
    % not DFT-normalised; multiplying by numT and halving converts them.
    idxVabcNode = salt.getRowIdxAbc_fromBusNum(busNum);
    vabc_SALT_1period = (1/2)*salt.getTimeDomainSignal_iDFT(idxVabcNode)*numT;
    vabc_SALT = repmat(vabc_SALT_1period(:,tsampIdxPeriodic_downsample), ...
        1,idx_periodic);

    figure, hold on
    % EMT solid and thick, SALT dashed in black on top: where SALT is
    % correct the dashed line disappears into the solid one.
    hTD = plot(tPlot,vabc_TD.','LineWidth',2.5);
    hSS = plot(tPlot,vabc_SALT(:,1:numTimepoints).','k--','LineWidth',1);
    legend([hTD(1) hTD(2) hTD(3) hSS(1)], ...
        'EMT phase a','EMT phase b','EMT phase c','SALT (steady state)', ...
        'Location','eastoutside','ItemHitFcn',@cb_legend)
    xlabel('Time [s]'), ylabel('Voltage [PU]')
    title("Bus " + num2str(busNum) + " - EMT time domain vs SALT steady state")
    grid on
    xlim([tPlot(1) tPlot(end)])

    exportgraphics(gcf,fullfile(resultsDir, ...
        sprintf('%s_bus%d.png',salt.caseName,busNum)),'Resolution',150);
end

end   % verifyWithEmtSim

%% Report the steady-state operating point
% The bus table is written to a file rather than the command window: it is
% one line per bus, which is 13659 lines on the largest case.
resultsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results');
if ~isfolder(resultsDir), mkdir(resultsDir), end
resultsFile = fullfile(resultsDir,[char(salt.caseName) '_results.txt']);

fid = fopen(resultsFile,'w');
assert(fid>0,'Could not open %s for writing',resultsFile)
closeResults = onCleanup(@() fclose(fid));

fprintf(fid,'================ SALT steady-state solution ================\n');
fprintf(fid,'Case             : %s\n',salt.caseName);
fprintf(fid,'System frequency : %.6f Hz  (w = %.6f rad/s)\n',fSys,wSys);
fprintf(fid,'Deviation from %g Hz : %+.6f Hz\n\n',salt.f0,fSys-salt.f0);
fprintf(fid,'%-8s %-16s %-16s\n','Bus','|V| [PU]','angle(V) [deg]');
fprintf(fid,'%-8s %-16s %-16s\n','---','--------','--------------');
for idxNode = 1:numel(salt.bus)
    idxV = salt.getRowIdxAbc_fromBusNum(salt.bus(idxNode));
    Vph_a = complex(salt.x_v(idxV(1)),salt.x_v(idxV(2)));
    fprintf(fid,'%-8d %-16.6f %-+16.4f\n',salt.bus(idxNode), ...
        abs(Vph_a),rad2deg(angle(Vph_a)));
end
fprintf(fid,'============================================================\n');
clear closeResults

fprintf('\n================ SALT steady-state solution ================\n');
fprintf('System frequency : %.6f Hz  (w = %.6f rad/s)\n',fSys,wSys);
fprintf('Deviation from %g Hz : %+.6f Hz\n',salt.f0,fSys-salt.f0);
fprintf('Bus voltages (%d buses) written to:\n  %s\n',numel(salt.bus),resultsFile);
fprintf('============================================================\n');
