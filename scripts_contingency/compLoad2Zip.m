function [Pzip_zip_PF,Qzip_zip_PF,Pzip_IM_PF,Qzip_IM_PF,Pzip_PF,Qzip_PF, ...
    compLoadZipBus_PF,IMpowerPercent_react] = compLoad2Zip(salt)
%COMPLOAD2ZIP_SS  Reduce each composite (IM+ZIP) load to ZIP coefficients.
%
%   [...] = compLoad2Zip_SS(salt) takes a CONVERGED SALT object and, for
%   every composite load, sweeps the terminal voltage magnitude around its
%   operating point, solves the isolated load device at each point, and
%   fits the resulting P(Vm) and Q(Vm) curves with a quadratic. The
%   quadratic is normalised into ZIP shares and blended with the static ZIP
%   part of the same load, giving the per-load ZIP rows a power flow with a
%   modified (per-bus) ZIP model can consume.
%
%   Outputs are [3 x numLoadComp] ZIP rows ordered (Z,I,P):
%     Pzip_zip_PF / Qzip_zip_PF   static-ZIP share of the load
%     Pzip_IM_PF  / Qzip_IM_PF    induction-motor share of the load
%     Pzip_PF     / Qzip_PF       the sum, which sums to 1 down each column
%     compLoadZipBus_PF           bus number of each column
%     IMpowerPercent_react        motor share of the load's Q (the active
%                                 share is a parameter, the reactive share
%                                 falls out of the initialisation)
%
%   Background: P of an induction motor is essentially constant in Vm (see
%   "Induction Motor Static Models for Power Flow and Voltage Stability
%   Studies", Carmona-Sanchez et al.), while Q varies, so Q is the term the
%   ZIP fit really has to work for. Only the "modified" power-flow method
%   is supported here: each load keeps its own ZIP coefficients, so no
%   averaging across loads is done.
%
%   The caller owns its own figures, so nothing is plotted and nothing is
%   closed here.

%% Set up

% The solved operating point is read straight off the SALT object; there is
% no case/device setup to redo.
nLoad = salt.numLoadComp;
w0 = salt.w0;

% ZIP coefficients fitted to the motor, per load, ordered (Z,I,P)
PZIP = zeros(3,nLoad);
QZIP = zeros(3,nLoad);
% Nominal PQ of the motor alone, and of the whole composite load
xPQ0_im = zeros(2,nLoad);
xPQ0 = zeros(2,nLoad);

% Local row/col pattern of one composite-load block. The kernel returns
% triplet VALUES only, so the dense per-device Jacobian this sweep needs is
% rebuilt from the same map SALTASSEMBLEY uses for the real matrix.
mapLC = salt.map.loadComp;
nE = salt.numLoadCompEqns;

%% Sweep each load's terminal voltage magnitude

for idxLoad = 1:nLoad
    % Parameters. loadCompIM_R/_L/_LTe are already stored as [4 x 4 x nLoad]
    % pages, so slicing one load keeps the page shape the kernel expects.
    IM_R = salt.loadCompIM_R(:,:,idxLoad);
    IM_L = salt.loadCompIM_L(:,:,idxLoad);
    IM_LTe = salt.loadCompIM_LTe(:,:,idxLoad);
    IM_Tm0 = salt.loadCompIM_Tm0(idxLoad);
    IM_m = salt.loadCompIM_m(idxLoad);
    % The static ZIP part is switched off: this sweep characterises the
    % MOTOR alone, and the ZIP part is added back analytically at the end.
    Zip_PQz = zeros(2,1,1);
    Zip_PQi = zeros(2,1,1);
    Zip_PQp = zeros(2,1,1);
    [~,idxBus] = ismember(salt.loadComp_bus(idxLoad),salt.bus);
    busToLoadCompV = salt.busToLoadComp_V(idxBus,idxLoad);
    busToLoadCompI = salt.busToLoadComp_I(idxBus,idxLoad);

    % Nominal states, read out of the converged x_v
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
        [real(iphA0);imag(iphA0);real(iphB0);imag(iphB0);real(iphC0);imag(iphC0)];

    % Perturbed terminal voltage keeps the operating-point phase and only
    % scales the magnitude
    phaseMult = complex(vphA0(1),vphA0(2))/abs(complex(vphA0(1),vphA0(2)));
    assert(abs(phaseMult*zipVm0 - complex(vphA0(1),vphA0(2)))<1e-6)
    calcVphA = @(Vm) [real(calc3Ph(phaseMult*Vm));imag(calc3Ph(phaseMult*Vm))];

    % (Sweep)
    VmSweep = zipVm0*linspace(0.6,1.4,100);

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

    % (Save PQ)
    PQ_IM = zeros(2,numel(VmSweep));

    % (Solve)
    totIters_NR = 5;
    relErrNR = 1e-4;
    absErrNR = 1e-8;
    for idx = 1:numel(VmSweep)
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
            IRHS(idxVphA) = calcVphA(VmSweep(idx));
            % System
            Y(idxWsys,idxWsys) = 1;
            IRHS(idxWsys) = w0;
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
        xPQ = (x(idxVs(1))+1j*x(idxVs(2)))*(x(idxIs(1))-1j*x(idxIs(2)));
        PQ_IM(:,idx) = [real(xPQ);imag(xPQ)];
    end

    % Save
    xPQ0_im_complex = (vs0(1)+1j*vs0(2))*(is0(1)-1j*is0(2));
    xPQ0_im(:,idxLoad) = [real(xPQ0_im_complex);imag(xPQ0_im_complex)];
    xPQ0(:,idxLoad) = xPQ0_im(:,idxLoad) + salt.x_v(idxRowLoadComp+idxPQ);

    % Fitting - P
    zipCoeffsP = polyfit(VmSweep,PQ_IM(1,:),2);
    PZIP(:,idxLoad) = zipCoeffsP.*[zipVm0^2 zipVm0 1]/xPQ0_im(1,idxLoad);
    % Fitting - Q
    zipCoeffsQ = polyfit(VmSweep,PQ_IM(2,:),2);
    QZIP(:,idxLoad) = zipCoeffsQ.*[zipVm0^2 zipVm0 1]/xPQ0_im(2,idxLoad);
end

%% Finalize

% Each load keeps its own fit.
Pzip_final = PZIP;
Qzip_final = QZIP;
% Check is close to summing up to 1
assert(all((abs(sum(Pzip_final)).'-1)<0.1))
assert(all((abs(sum(Qzip_final)).'-1)<0.1))
% Make sum up to 1
% NOTE: The normalised ZIP is not identical to the raw quadratic fit. The
% normalisation forces P and Q to be exactly right at nominal voltage,
% whereas the unconstrained fit only minimises L2 error over the sweep.
Pzip_final = Pzip_final./sum(Pzip_final);
Qzip_final = Qzip_final./sum(Qzip_final);

% Static ZIP shares of the same loads, taken from the SALT object rather
% than from caller arguments. loadCompZip_PQabc is [2 x numLoadComp x 3]:
% dim 1 is the P/Q row, dim 2 is the load, dim 3 is (Z,I,P) - see
% compLoadValInitCalculations_assumeW0:
% >> PQz = PQ_Zip.*(obj.loadCompZip_PQabc(:,:,1)./(Vm.^2));
% >> PQi = PQ_Zip.*(obj.loadCompZip_PQabc(:,:,2)./Vm);
% >> PQp = PQ_Zip.*obj.loadCompZip_PQabc(:,:,3);
% PERMUTE puts (Z,I,P) on the rows and the load on the columns, matching
% the [3 x numLoadComp] orientation of the fitted motor shares.
PzipStatic = permute(salt.loadCompZip_PQabc(1,:,:),[3 2 1]);
QzipStatic = permute(salt.loadCompZip_PQabc(2,:,:),[3 2 1]);

% Motor share of the load's active power. This is per load in SALT, so it
% is used as a [1 x numLoadComp] row and broadcast down the three ZIP rows.
IMpowerPercent = salt.loadCompIM_powerPercent(:).';

% Write output
% NOTE: These parameters do NOT need to worry about bus/IM PU because these
% are ratios, not coefficients which will be calculated later!!
assert(all(abs(IMpowerPercent-xPQ0_im(1,:)./xPQ0(1,:))<1e-6))
IMpowerPercent_react = xPQ0_im(2,:)./xPQ0(2,:);
%
Pzip_zip_PF = (1-IMpowerPercent).*PzipStatic;
Qzip_zip_PF = (1-IMpowerPercent_react).*QzipStatic;
%
Pzip_IM_PF = IMpowerPercent.*Pzip_final;
Qzip_IM_PF = IMpowerPercent_react.*Qzip_final;
%
Pzip_PF = Pzip_zip_PF + Pzip_IM_PF;
Qzip_PF = Qzip_zip_PF + Qzip_IM_PF;
%
assert(all(abs(sum(Pzip_PF)-1)<1e-6))
assert(all(abs(sum(Qzip_PF)-1)<1e-6))

compLoadZipBus_PF = salt.loadComp_bus;

end
