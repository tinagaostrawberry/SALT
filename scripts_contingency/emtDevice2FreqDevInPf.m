function [dP_dw_gen,dP_dw_ibrGFM,dP_dw_load,dQ_dw_load] = ...
    emtDevice2FreqDevInPf(salt)
%EMTDEVICE2FREQDEVINPF_SS  Frequency droop of every device, for the power flow.
%
%   [...] = emtDevice2FreqDevInPf_SS(salt) takes a CONVERGED SALT object and
%   returns, per device, the sensitivity of injected power to a deviation
%   in system frequency, in the form a frequency-aware power flow consumes.
%
%     dP_dw_gen     [numGenrou x 2]   [bus, dP/dw] of each GENROU, from its
%                                     governor droop and damping
%     dP_dw_ibrGFM  [numIbrGFM x 2]   [bus, dP/dw] of each GFM IBR, from its
%                                     active-power droop gain Mp
%     dP_dw_load    [numLoadComp x 2] [bus, dP/dw] of each composite load
%     dQ_dw_load    [numLoadComp x 2] [bus, dQ/dw] of each composite load
%
%   Machine and IBR droops are analytic. The composite loads have no closed
%   form, so each one is solved in isolation over a frequency sweep and the
%   result is fitted with the linear model the power flow assumes,
%     P = P0*(1 + dP_dw*deltaW),  Q = Q0*(1 + dQ_dw*deltaW).
%
%   The caller owns its own figures, so nothing is plotted and nothing is
%   closed here.

%% Generators

% Pmech = Pt + Ra*It^2
% (Pref-deltaWr_PU)/R - deltaWr_PU*Dt = Pmech
% =>
% Pt = Pref/R - deltaWr_PU*(1/R + Dt) ...ignoring Ra

% In SALT the droop parameters live on the TGOV1D devices, not on the
% machines, and tgov1d_genrouIdx says which GENROU each governor drives.
% Build a per-GENROU droop vector instead of assuming governor k belongs to
% machine k, which is only true when the device file happens to list them
% in the same order.
govGenrouIdx = salt.tgov1d_genrouIdx(:);
govPerGenrou = accumarray(govGenrouIdx,1,[salt.numGenrou,1]);
if any(govPerGenrou~=1)
    error(['emtDevice2FreqDevInPf_SS: every GENROU needs exactly one ' ...
        'TGOV1D governor to define its droop, but GENROU(s) [%s] have ' ...
        'a different count'],num2str(find(govPerGenrou~=1).'))
end
gov_R = zeros(salt.numGenrou,1);
gov_Dt = zeros(salt.numGenrou,1);
gov_R(govGenrouIdx) = salt.tgov1d_R(:);
gov_Dt(govGenrouIdx) = salt.tgov1d_Dt(:);

% Force both index lists to COLUMNS. The per-unit conversions are
% elementwise, so a column bus index against a row device index would
% broadcast into an n-by-n matrix instead of the n-by-1 list wanted here.
[~,bus_idx] = ismember(salt.genrou_bus,salt.bus);
bus_idx = bus_idx(:);
gen_idx = (1:salt.numGenrou).';
genToBusV = salt.genrouToBus_V(bus_idx,gen_idx);
genToBusI = salt.genrouToBus_I(bus_idx,gen_idx);

dP_dw_gen = (-(1./gov_R + gov_Dt)/salt.w0) .* (genToBusV.*genToBusI);
dP_dw_gen = [salt.genrou_bus(:),dP_dw_gen(:)];

%% GF IBRs

% wr = w0 + Mp*(pset - pavg)
% =>
% p = pset + deltaWr(-1/Mp)

[~,bus_idx] = ismember(salt.ibrGFM_bus,salt.bus);
bus_idx = bus_idx(:);                       % columns, as above
ibrGFM_idx = (1:salt.numIbrGFM).';
ibrGFMToBusV = salt.ibrGFMToBus_V(bus_idx,ibrGFM_idx);
ibrGFMToBusI = salt.ibrGFMToBus_I(bus_idx,ibrGFM_idx);

dP_dw_ibrGFM = (-1./salt.ibrGFM_Mp(:)) .* (ibrGFMToBusV.*ibrGFMToBusI);
dP_dw_ibrGFM = [salt.ibrGFM_bus(:),dP_dw_ibrGFM(:)];

%% Loads

% Sweep of frequency
f0 = salt.w0/(2*pi);
w_sweep = 2*pi*(f0 + ((-2):0.1:(2)));
% Save P,Q values w.r.t frequency
PQ_sweep = zeros(2,numel(w_sweep));
% Save dP,dQ w.r.t. frequency
dP_dw_load = zeros(salt.numLoadComp,1);
dQ_dw_load = zeros(salt.numLoadComp,1);

% Local row/col pattern of one composite-load block. The kernel returns
% triplet VALUES only, so the dense per-device Jacobian this sweep needs is
% rebuilt from the same map SALTASSEMBLEY uses for the real matrix.
mapLC = salt.map.loadComp;
nE = salt.numLoadCompEqns;

% Calculate SS IM PQ w.r.t frequency
for idxLoad = 1:salt.numLoadComp

    % Load constants for this load. loadCompIM_R/_L/_LTe are already stored
    % as [4 x 4 x numLoadComp] pages, so slicing one load keeps the page
    % shape the kernel expects.
    IM_R = salt.loadCompIM_R(:,:,idxLoad);
    IM_L = salt.loadCompIM_L(:,:,idxLoad);
    IM_LTe = salt.loadCompIM_LTe(:,:,idxLoad);
    IM_Tm0 = salt.loadCompIM_Tm0(idxLoad);
    IM_m = salt.loadCompIM_m(idxLoad);
    % The static ZIP part is switched off: this sweep characterises the
    % MOTOR alone, which is the only frequency-dependent part of the load.
    Zip_PQz = zeros(2,1,1);
    Zip_PQi = zeros(2,1,1);
    Zip_PQp = zeros(2,1,1);
    w0 = salt.w0;
    [~,idxBus] = ismember(salt.loadComp_bus(idxLoad),salt.bus);
    busToLoadCompV = salt.busToLoadComp_V(idxBus,idxLoad);
    busToLoadCompI = salt.busToLoadComp_I(idxBus,idxLoad);

    % Get initial/nominal values, read out of the converged x_v
    idxRowLoadComp = salt.getRowIdx_fromLoadCompIndices(idxLoad)-1;
    [idxVs,idxIphABC,idxWSlip,idxThetaSys,idxIs,idxIr,idxWr, ...
        idxPQ,idxVm,idxIzip] = getLoadCompIndices_SS;
    idxVphAbc = salt.getRowIdxAbc_fromBusNum(salt.loadComp_bus(idxLoad));
    vs0 = salt.x_v(idxRowLoadComp+idxVs);
    wslip0 = salt.x_v(idxRowLoadComp+idxWSlip);
    thetaSys0 = salt.x_v(idxRowLoadComp+idxThetaSys);
    is0 = salt.x_v(idxRowLoadComp+idxIs);
    ir0 = salt.x_v(idxRowLoadComp+idxIr);
    wr0 = salt.x_v(idxRowLoadComp+idxWr);
    vphA0 = salt.x_v(idxVphAbc(1:2));
    zipPQ0 = zeros(2,1);
    zipVm0 = salt.x_v(idxRowLoadComp+idxVm);
    zipI = zeros(2,1);
    [iphA0,iphB0,iphC0] = calc3Ph(is0(1)+1j*is0(2));
    iphAbc0 = (1/busToLoadCompI)* ...
        [real(iphA0);imag(iphA0); ...
        real(iphB0);imag(iphB0); ...
        real(iphC0);imag(iphC0)];
    %
    PQ0_IM_complex = (vs0(1)+1j*vs0(2))*(is0(1)-1j*is0(2));
    P0_IM = real(PQ0_IM_complex);
    Q0_IM = imag(PQ0_IM_complex);

    % Calculate P,Q w.r.t. frequency
    % (Initialize) the device block plus the two terminal-voltage rows and
    % the system-frequency row that pin the sweep
    matSize = nE+2+1;
    x0 = zeros(matSize,1);
    idxVphA = nE+(1:2);
    idxWsys = nE+3;
    x0( ...
        [idxVs,idxIphABC,idxWSlip,idxThetaSys,idxIs,idxIr,idxWr,idxPQ,idxVm,idxIzip, ...
        idxVphA, ...
        idxWsys]) = ...
        [vs0;iphAbc0;wslip0;thetaSys0;is0;ir0;wr0;zipPQ0;zipVm0;zipI; ...
        vphA0; ...
        w0];
    % (Solve)
    totIters_NR = 5;
    relErrNR = 1e-4;
    absErrNR = 1e-8;
    for idxSweep = 1:numel(w_sweep)
        xprev = x0;
        x = x0;
        for idxNR = 1:totIters_NR
            % Get matrix stamps. Every state/parameter goes in as an
            % [nState x 1 x n] page with n = 1, because the kernel is
            % vectorised over loads and infers n from size(wr,3).
            [Jvals,FloadComp, ...
                dFParkV_dVphA,dFslip_dwSys,dFs_dwSys,dFswing_dwSys] = ...
                stampOfCompositeLoadMatrices_SS( ...
                ...% States
                reshape(x(idxVs),2,1,1),reshape(x(idxIphABC),6,1,1), ...
                reshape(x(idxWSlip),1,1,1),reshape(x(idxThetaSys),1,1,1), ...
                reshape(x(idxIs),2,1,1),reshape(x(idxIr),2,1,1), ...
                reshape(x(idxWr),1,1,1),reshape(x(idxVphA),2,1,1), ...
                x(idxWsys), ...
                reshape(x(idxPQ),2,1,1),reshape(x(idxVm),1,1,1), ...
                reshape(x(idxIzip),2,1,1), ...
                ... % Parameters
                reshape(IM_R,4,4,1),reshape(IM_L,4,4,1), ...
                reshape(IM_LTe,4,4,1),reshape(IM_Tm0,1,1,1), ...
                reshape(IM_m,1,1,1),Zip_PQz,Zip_PQi,Zip_PQp, ...
                ... % Consts
                w0,busToLoadCompI,busToLoadCompV);
            % Rebuild the dense block Jacobian from the triplet values.
            % ACCUMARRAY (not sparse assignment) because several stamps
            % land on the same (row,col) and must add, exactly as
            % SALTASSEMBLEY accumulates them in the real system matrix.
            if numel(mapLC.row)~=size(Jvals,1)
                error(['Composite-load stamp map has %d triplets but the ' ...
                    'kernel returned %d; map and kernel are out of sync'], ...
                    numel(mapLC.row),size(Jvals,1))
            end
            JloadComp = accumarray([mapLC.row(:) mapLC.col(:)],Jvals(:),[nE nE]);
            % The Park coupling comes back flattened as [4 x n]
            dParkV = reshape(dFParkV_dVphA,2,2);
            % Create matrix
            Y = zeros(matSize);
            IRHS = zeros(matSize,1);
            % Stamp in nonlinear portion
            % J(x[k])*x[k+1] = J(x[k])*x[k]-F(x[k])
            Y(1:nE,1:nE) = JloadComp;
            IRHS(1:nE) = JloadComp*x(1:nE) - FloadComp;
            % Park
            Y(idxVs,idxVphA) = dParkV;
            IRHS(idxVs) = IRHS(idxVs) + dParkV*x(idxVphA);
            % Slip
            Y(idxWSlip,idxWsys) = dFslip_dwSys;
            IRHS(idxWSlip) = IRHS(idxWSlip) + dFslip_dwSys*x(idxWsys);
            % Stator
            Y(idxIs,idxWsys) = dFs_dwSys;
            IRHS(idxIs) = IRHS(idxIs) + dFs_dwSys*x(idxWsys);
            % Swing
            Y(idxWr,idxWsys) = dFswing_dwSys;
            IRHS(idxWr) = IRHS(idxWr) + dFswing_dwSys*x(idxWsys);
            % Voltage source
            Y(idxVphA,idxVphA) = eye(2);
            IRHS(idxVphA) = vphA0;
            % System
            Y(idxWsys,idxWsys) = 1;
            IRHS(idxWsys) = w_sweep(idxSweep);
            % Solve
            x = Y \ IRHS;
            % Check convergence
            converge = isErrSmall(x,xprev,relErrNR,absErrNR);
            if converge
                break
            else
                if (idxNR<totIters_NR)
                    xprev = x;
                end
            end
        end
        if ~converge
            error('Did not converge')
        end
        % Save
        vs = x(idxVs);
        is = x(idxIs);
        PQ_complex = (vs(1)+1j*vs(2))*(is(1)-1j*is(2));
        PQ_sweep(:,idxSweep) = [real(PQ_complex);imag(PQ_complex)];
    end

    % Fitting
    % P = (Pz*Vm^2 + Pi*Vm + Pp) * (1 + dP*deltaW)
    % Q = (Qz*Vm^2 + Qi*Vm + Qp) * (1 + dQ*deltaW)
    dP_IM = (w_sweep.'-w0) \ ((PQ_sweep(1,:).')/P0_IM - 1);
    dQ_IM = (w_sweep.'-w0) \ ((PQ_sweep(2,:).')/Q0_IM - 1);

    % Saving
    % NOTE: Do NOT need to normalize PU because equation is:
    % >> P = P0*(1+dP_dw*deltaW)
    % ... so any scaling that happens to P would happen to P0 as a
    % multiplier
    dP_dw_load(idxLoad) = dP_IM;
    dQ_dw_load(idxLoad) = dQ_IM;
end

dP_dw_load = [salt.loadComp_bus(:),dP_dw_load];
dQ_dw_load = [salt.loadComp_bus(:),dQ_dw_load];

% Without composite loads there is nothing frequency-dependent on the load
% side, so report the plain constant-impedance loads with zero droop. This
% replaces the legacy include_compositeLoads flag, which the caller no
% longer carries.
if salt.numLoadComp == 0
    dP_dw_load = [salt.load_bus(:),zeros(numel(salt.load_bus),1)];
    dQ_dw_load = [salt.load_bus(:),zeros(numel(salt.load_bus),1)];
end

end
