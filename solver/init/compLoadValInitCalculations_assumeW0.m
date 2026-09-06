function [vdq_s,loadComp_IComplexPhABC,wslip,thetaSys,wr,idq_s,idq_r, ...
    PQ_Zip,Vm,loadComp_IComplexPhA_Zip,T0,PQz,PQi,PQp] = ...
    compLoadValInitCalculations_assumeW0(obj,bus_VComplex)
%COMPLOADVALINITCALCULATIONS_ASSUMEW0  Composite-load operating point.
%
%   A composite load is one induction motor and one ZIP load sharing a bus.
%   The power flow gives the bus its terminal voltage and a single dispatch
%   P + jQ; this routine divides that between the two and finds the motor's
%   internal state, assuming nominal frequency.
%
%   1. PER-UNIT BASE, AND THE ONLY PLACE THE POWER FACTOR IS USED. The
%      motor's resistances and reactances are given per unit on its own
%      rating, so a rating has to be chosen consistent with them. It is
%      taken as the motor's apparent power at full load,
%
%          IM_BaseMVA = P_bus * powerPercent / PF
%
%      where powerPercent is the motor's share of the bus real power and PF
%      is its assumed power factor. PF converts the motor's real power into
%      an apparent power and so sets the base; it does NOT set the motor's
%      reactive power, and nothing downstream reads it again. Bus voltage,
%      current and dispatch are converted onto this base first.
%
%   2. MOTOR REAL POWER IS ASSIGNED, REACTIVE POWER IS NOT. The motor takes
%      powerPercent of the bus real power. Its reactive power is whatever
%      the machine equations produce at the slip that carries that real
%      power, which is generally not what the power factor above would
%      suggest. A first guess Q = sqrt(S^2 - P^2) at S = 1 pu (fully loaded)
%      only seeds the Newton iteration below.
%
%   3. NEWTON SOLVE. The stator and rotor dq currents, the slip, the rotor
%      speed, the mechanical torque and the resulting reactive power are
%      solved together from the motor's steady-state equations.
%
%   4. THE ZIP TAKES THE REMAINDER. Real power is split by fraction,
%      P_zip = (1 - powerPercent)*P_bus, and reactive power by SUBTRACTION,
%      Q_zip = Q_bus - Q_IM. Subtraction rather than a fraction is what lets
%      a bus with zero (or negative) reactive dispatch initialise: the ZIP
%      simply supplies -Q_IM. The remainder is then distributed over the
%      constant-impedance, constant-current and constant-power coefficients
%      at the initial voltage. See the EPRI note in the ZIP section below
%      for why the reactive split is left to come out of the solve rather
%      than being forced to a specified power factor.
%
%   5. CHECKS. The motor and ZIP currents must add up to the bus current,
%      and their apparent powers to the bus dispatch, at the load base and
%      back at the bus base.

[~,busIdx_loadComp] = ismembertol(obj.loadComp_bus,obj.bus);

% Overall load
% Process: Solve in bus PU, convert to load PU
% 1. Calc voltage and current (line-to-line (LL), RMS, bus PU):
loadComp_VComplexPhA_busPU = bus_VComplex(busIdx_loadComp);
loadComp_SComplex_busPU = (obj.loadComp_MW + 1j*obj.loadComp_MVar)/ ...
    obj.baseMVA;
loadComp_IComplexPhA_busPU = ((loadComp_SComplex_busPU')./ ...
    (loadComp_VComplexPhA_busPU')).';
% 2. Convert from bus or raw (MW,MVAR) to load PU
loadComp_IComplexPhA = loadComp_IComplexPhA_busPU.* ...
    obj.busToLoadComp_I(busIdx_loadComp,1:obj.numLoadComp);
loadComp_VComplexPhA = loadComp_VComplexPhA_busPU.* ...
    obj.busToLoadComp_V(busIdx_loadComp,1:obj.numLoadComp);
loadComp_SAbs = abs((obj.loadComp_MW + 1j*obj.loadComp_MVar)./ ...
    obj.loadComp_MVA_base);
assert(all(abs(loadComp_SAbs-abs(loadComp_VComplexPhA.*conj(loadComp_IComplexPhA))) ...
    <1e-6))

%% Initialize induction motor

% Calculate IM's P in PU - is fraction of entire load
P_IM = obj.loadCompIM_powerPercent.*obj.loadComp_MW ./ ...
    obj.loadComp_MVA_base;
% Estimate MVA as fully loaded (in PU, should equal to 1)
loadComp_SAbs_IM_est = ones(obj.numLoadComp,1);
if any(loadComp_SAbs_IM_est<P_IM)
    % However, if due to load perturbation (which doesn't change base MVA
    % but does change the actual operating MVA), it is more than fully
    % loaded
    idxIncreaseLoad = loadComp_SAbs_IM_est<P_IM;
    loadComp_SAbs_IM_est(idxIncreaseLoad) = ...
        (P_IM(idxIncreaseLoad).^2 + 0.01).^(1/2);
end
% Estimate Q and the corresponding current
Q_IM_est = sqrt(loadComp_SAbs_IM_est.^2-P_IM.^2);
loadComp_IComplexPhA_IM_est = ...
    (complex(P_IM,Q_IM_est)' ./ loadComp_VComplexPhA').';


% Use solution based on PQ estimate, which will be inaccurate (violate q
% rotor eqn, specifically), as initial guess
% (Calc stator)
vdq_s0 = [real(loadComp_VComplexPhA).';imag(loadComp_VComplexPhA).'];
idq_s0 = [real(loadComp_IComplexPhA_IM_est).'; ...
    imag(loadComp_IComplexPhA_IM_est).'];
% (Est rotor and iabc)
freqMat = [0 -1 0 0; 1 0 0 0]; % (1/w0)[0 -wsys 0 0;wsys 0 0 0];
idq_r0 = zeros(2,obj.numLoadComp);
iabc0 = zeros(3,obj.numLoadComp);
for idx = 1:obj.numLoadComp
    Y = freqMat*obj.loadCompIM_L(:,3:4,idx);
    I = vdq_s0(:,idx) - ...
        (obj.loadCompIM_R(1:2,1:2,idx)+freqMat*obj.loadCompIM_L(:,1:2,idx)) ...
        *idq_s0(:,idx);
    idq_r0(:,idx) = Y\I;
    %
    [iphA0,iphB0,iphC0] = calc3Ph(loadComp_IComplexPhA_busPU(idx));
    iabc0(:,idx) = real([iphA0;iphB0;iphC0]);
end
% (Calc wslip)
Rr = squeeze(obj.loadCompIM_R(3,3,:)).';
Lm = squeeze(obj.loadCompIM_L(1,3,:)).';
Lr = squeeze(obj.loadCompIM_L(3,3,:)).';
wslip0 = obj.w0 * ( Rr.*idq_r0(1,:) ) ./ ...
    (Lm.*idq_s0(2,:) + Lr.*idq_r0(2,:));
% (Calc wr)
wr0 = obj.w0 - wslip0;
% (Calc torque constant - make efficient without for loop?)
idq_sr0 = reshape([idq_s0;idq_r0],4,1,obj.numLoadComp);
T00 = ...
    squeeze( ...
    pagemtimes(pagemtimes(pagetranspose(idq_sr0),obj.loadCompIM_LTe),idq_sr0) ...
    ).' ...
    ./ (wr0/obj.w0).^(obj.loadCompIM_m.');

% Solve for SS solution
% Initialize solution
vdq_s = zeros(2,obj.numLoadComp);
wslip = zeros(1,obj.numLoadComp);
thetaSys = zeros(1,obj.numLoadComp);
idq_s = zeros(2,obj.numLoadComp);
idq_r = zeros(2,obj.numLoadComp);
wr = zeros(1,obj.numLoadComp);
T0 = zeros(1,obj.numLoadComp);
PQ_IM = [P_IM.';zeros(1,obj.numLoadComp)];
% Indices
idxVs_ss = 1:2;
idxWSlip_ss = 3;
idxIs_ss = 4:5;
idxIr_ss = 6:7;
idxWr_ss = 8;
idxT0_ss = 9;
idxQ_ss = 10;
% Iteratively solve
for idx = 1:obj.numLoadComp
    % Pre-allocate constants
    R = obj.loadCompIM_R(:,:,idx);
    L_PU = obj.loadCompIM_L(:,:,idx)/obj.w0;
    L_Te = obj.loadCompIM_LTe(:,:,idx);
    m = obj.loadCompIM_m(idx);
    P = obj.loadCompIM_powerPercent(idx)*obj.loadComp_MW(idx) / ...
        obj.loadComp_MVA_base(idx);

    % Pre-allocate equations
    % (Park V)
    FParkV = @(x) -x(idxVs_ss) + [real(loadComp_VComplexPhA(idx)); ...
        imag(loadComp_VComplexPhA(idx))];
    dFParkV_dVs = @(x) -eye(2);
    % (Slip frequency)
    Fslip = @(x) -x(idxWSlip_ss) + obj.w0 - x(idxWr_ss);
    dFslip_dwSlip_dwr = @(~) -ones(1,2);
    % (Stator/rotor)
    wMat = @(x) [0 -obj.w0 0 0;obj.w0 0 0 0;0 0 0 -x(idxWSlip_ss); ...
        0 0 x(idxWSlip_ss) 0];
    phisr = @(x,eqnIdx) L_PU(eqnIdx,:)*x([idxIs_ss,idxIr_ss]);
    Fsr = @(x) -[x(idxVs_ss);zeros(2,1)]+R*x([idxIs_ss,idxIr_ss])+ ...
        wMat(x)*phisr(x,1:4);
    dFsr_dVs_dIsr_dWSlip = @(x) [-[eye(2);zeros(2)],R+wMat(x)*L_PU, ...
        [0;0;-phisr(x,4);phisr(x,3)]];
    % (Swing)
    Fswing = @(x) ...
        (x([idxIs_ss,idxIr_ss]).')*L_Te*x([idxIs_ss,idxIr_ss]) ...
        - x(idxT0_ss)*(x(idxWr_ss)/obj.w0)^m;
    dFswing_dWr_dIsr_dT0 = @(x) [-x(idxT0_ss)*(m*x(idxWr_ss)^(m-1))/(obj.w0^m), ...
        x([idxIs_ss,idxIr_ss]).'*(L_Te+L_Te.').',-(x(idxWr_ss)/obj.w0)^m];
    % (P)
    FP = @(x) -P + x(idxVs_ss).'*x(idxIs_ss);
    dFP_dVs_dIs = @(x) [x(idxIs_ss).',x(idxVs_ss).'];
    % (Q)
    FQ = @(x) -x(idxQ_ss) + x(idxIs_ss).'*[0 1;-1 0]*x(idxVs_ss);
    dFQ_dQ_dVs_dIs = @(x) [-1,x(idxIs_ss).'*[0 1;-1 0], ...
        x(idxVs_ss).'*[0 1;-1 0].'];

    % Initialize
    x = zeros(10,1);
    x([idxVs_ss,idxWSlip_ss,idxIs_ss,idxIr_ss,idxWr_ss,idxT0_ss,idxQ_ss]) = ...
        [vdq_s0(:,idx);wslip0(idx);idq_s0(:,idx);idq_r0(:,idx);wr0(idx); ...
        T00(idx);idq_s0(:,idx).'*[0 1;-1 0]*vdq_s0(:,idx)];

    % Iteratively solve
    totIters_NR = 500;
    relErrNR = 1e-4;
    absErrNR = 1e-8;
    for idxNR = 1:totIters_NR
        % Create matrix
        Y = zeros(10);
        I = zeros(10,1);
        Y(1:2,idxVs_ss) = dFParkV_dVs(x);
        I(1:2) = dFParkV_dVs(x)*x(idxVs_ss) - FParkV(x);
        Y(3,[idxWSlip_ss,idxWr_ss]) = dFslip_dwSlip_dwr(x);
        I(3) = dFslip_dwSlip_dwr(x)*x([idxWSlip_ss,idxWr_ss]) - Fslip(x);
        Y(4:7,[idxVs_ss,idxIs_ss,idxIr_ss,idxWSlip_ss]) = dFsr_dVs_dIsr_dWSlip(x);
        I(4:7) = dFsr_dVs_dIsr_dWSlip(x)* ...
            x([idxVs_ss,idxIs_ss,idxIr_ss,idxWSlip_ss]) - Fsr(x);
        Y(8,[idxWr_ss,idxIs_ss,idxIr_ss,idxT0_ss]) = dFswing_dWr_dIsr_dT0(x);
        I(8) = dFswing_dWr_dIsr_dT0(x)*x([idxWr_ss,idxIs_ss,idxIr_ss,idxT0_ss]) ...
            - Fswing(x);
        Y(9,[idxVs_ss,idxIs_ss]) = dFP_dVs_dIs(x);
        I(9) = dFP_dVs_dIs(x)*x([idxVs_ss,idxIs_ss]) - FP(x);
        Y(10,[idxQ_ss,idxVs_ss,idxIs_ss]) = dFQ_dQ_dVs_dIs(x);
        I(10) = dFQ_dQ_dVs_dIs(x)*x([idxQ_ss,idxVs_ss,idxIs_ss]) - FQ(x);
        % Solve update
        x_update = Y\I;
        %x_update = stepLimiting(x_update,x,[0 0.5]);
        % Converge?
        converge = isErrSmall(x_update,x,relErrNR,absErrNR);
        if converge
            break
        else
            if (idxNR<totIters_NR)
                x = x_update;
            end
        end
    end
    if ~converge
        error('Initialization of IM in composite load failed!')
    end

    % Save
    vdq_s(:,idx) = x(idxVs_ss);
    wslip(idx) = x(idxWSlip_ss);
    thetaSys(idx) = 0;
    idq_s(:,idx) = x(idxIs_ss);
    idq_r(:,idx) = x(idxIr_ss);
    wr(idx) = x(idxWr_ss);
    T0(idx) = x(idxT0_ss);
    PQ_IM(2,idx) = x(idxQ_ss);

end

%% ZIP
% Calculate PQ
% According to EPRI's "Technical Reference on the Composite Load Model",
% the reactive power (Q) is ...
% "...computed based on the user specified power factor. It is  possible
% that the total reactive power computed during initialization does not
% match with the load flow solution. In this case, a computed additional
% shunt capacitance, Bfdr, is added to the composite load model to account
% for this mismatch."
% However, HERE, we do not specify the power factor to force the value of
% the ZIP, but just simply allow the remaining Q to be absorbed by ZIP...
PQ_Zip = [ ...
    ... % P
    ((1-obj.loadCompIM_powerPercent).* ...
    obj.loadComp_MW./obj.loadComp_MVA_base).'; ...
    ... % Q
    (obj.loadComp_MVar./obj.loadComp_MVA_base).'-PQ_IM(2,:) ...
    ];

% Calc voltage magnitude
Vm = abs(loadComp_VComplexPhA).';

% Calculate constants
PQz = PQ_Zip.*(obj.loadCompZip_PQabc(:,:,1)./(Vm.^2));
PQi = PQ_Zip.*(obj.loadCompZip_PQabc(:,:,2)./Vm);
PQp = PQ_Zip.*obj.loadCompZip_PQabc(:,:,3);

% Solve for current in dq frame
loadComp_IComplexPhA_IM = ...
    (complex(PQ_IM(1,:),-PQ_IM(2,:)) ./ ...
    loadComp_VComplexPhA').';
loadComp_IComplexPhA_Zip = ...
    (complex(PQ_Zip(1,:),-PQ_Zip(2,:)) ./ ...
    loadComp_VComplexPhA').';

%% Total current

% Place in network
loadComp_IComplexPhABC = zeros(3,obj.numLoadComp);
for idxLoad = 1:obj.numLoadComp
    [iphA,iphB,iphC] = calc3Ph(loadComp_IComplexPhA_IM(idxLoad) + ...
        loadComp_IComplexPhA_Zip(idxLoad));
    [~,idxBus] = ismember(obj.loadComp_bus(idxLoad),obj.bus);
    iphABC = [iphA;iphB;iphC] * ...
        obj.loadCompToBus_I(idxBus,idxLoad);
    loadComp_IComplexPhABC(:,idxLoad) = iphABC;
end


%% IM-ZIP unity check

% (Current)
totCurr = loadComp_IComplexPhA - loadComp_IComplexPhA_IM ...
    - loadComp_IComplexPhA_Zip;
assert(all(abs(totCurr)<1e-6))
% (Power, IM)
S_IM = complex(PQ_IM(1,:),PQ_IM(2,:)).';
S_calc_IM = loadComp_VComplexPhA .* ...
    conj(loadComp_IComplexPhA_IM);
assert(all(abs(S_IM-S_calc_IM)<1e-6))
% (Power, ZIP)
S_ZIP = complex(PQ_Zip(1,:),PQ_Zip(2,:)).';
S_calc_ZIP = loadComp_VComplexPhA .* ...
    conj(loadComp_IComplexPhA_Zip);
assert(all(abs(S_ZIP-S_calc_ZIP)<1e-6))
% (Power, total)
SComplex = loadComp_SComplex_busPU*obj.baseMVA ./ obj.loadComp_MVA_base;
assert(all(abs(SComplex-S_IM-S_ZIP)<1e-6))
% (Power, total; convert to bus PU, solve in bus PU)
loadComp_IComplexPhA_IM_busPU = ...
    loadComp_IComplexPhA_IM .* ...
    obj.loadCompToBus_I(busIdx_loadComp,1:obj.numLoadComp);
loadComp_IComplexPhA_Zip_busPU = ...
    loadComp_IComplexPhA_Zip .* ...
    obj.loadCompToBus_I(busIdx_loadComp,1:obj.numLoadComp);
S_calc_IM_busPU = loadComp_VComplexPhA_busPU .* ...
    conj(loadComp_IComplexPhA_IM_busPU);
S_calc_ZIP_busPU = loadComp_VComplexPhA_busPU .* ...
    conj(loadComp_IComplexPhA_Zip_busPU);
assert(all(abs(loadComp_SComplex_busPU-S_calc_IM_busPU ...
    -S_calc_ZIP_busPU)<1e-6))

