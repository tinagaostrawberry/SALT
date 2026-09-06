function [obj,T0,PQz,PQi,PQp] = loadCompValInitializations(obj, ...
    bus_VComplex)

[~,busIdx_loadComp] = ismembertol(obj.loadComp_bus,obj.bus);

% Overall load
% Process: Solve in bus PU, convert to load PU
% 1. Calc voltage and current (line-to-line (LL), RMS, bus PU):
loadComp_VComplexPhA_busPU = bus_VComplex(busIdx_loadComp);
loadComp_SComplex_busPU = (obj.loadComp_MW + 1j*obj.loadComp_MVar)/ ...
    obj.baseMVA;
loadComp_IComplexPhA_busPU = ((loadComp_SComplex_busPU')./ ...
    (loadComp_VComplexPhA_busPU')).';
% 2. Convert from bus to load PU
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
% Estimate Q and the corresponding current
Q_IM_est = sqrt(loadComp_SAbs.^2-P_IM.^2);
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
%assert(all(wslip0/obj.w0<0.05))
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
    totIters_NR = 500; % Step limiting
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
        %x_update = stepLimiting(x_update,x,[0 1]);
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

% Check slip frequency is within 95% of nominal
%assert(all(wslip0/obj.w0<0.05))
% Check power factor is around 0.85 (within 0.1)
%PF_IM = PQ_IM(1,:) ./ sqrt(PQ_IM(1,:).^2 + PQ_IM(2,:).^2);
%assert(all(abs(PF_IM-0.85)<0.1))

%% ZIP
% Calculate PQ
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

%% Total current

% Solve for current in dq frame
loadComp_IComplexPhA_IM = ...
    (complex(PQ_IM(1,:),-PQ_IM(2,:)) ./ ...
    loadComp_VComplexPhA').';
loadComp_IComplexPhA_Zip = ...
    (complex(PQ_Zip(1,:),-PQ_Zip(2,:)) ./ ...
    loadComp_VComplexPhA').';

% Place in network
iabc = zeros(3,obj.numLoadComp);
for idxLoad = 1:obj.numLoadComp
    [iphA,iphB,iphC] = calc3Ph(loadComp_IComplexPhA_IM(idxLoad) + ...
        loadComp_IComplexPhA_Zip(idxLoad));
    [~,idxBus] = ismember(obj.loadComp_bus(idxLoad),obj.bus);
    iabc(:,idxLoad) = real([iphA;iphB;iphC]) * ...
        obj.loadCompToBus_I(idxBus,idxLoad);
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


%% Set into x_v
% Final initialize
rowIdxLoadComp = obj.getRowIdx_fromLoadCompIndices(1:obj.numLoadComp)-1;
[idxVs,idxIabc,idxWSlip,idxThetaSys,idxIs,idxIr,idxWr, ...
    idxPQ,idxVmPreFilt,idxVm,idxIzip] = getLoadCompIndicesTD;
%
obj.x_v(rowIdxLoadComp+idxVs(1)) = vdq_s(1,:);
obj.x_v(rowIdxLoadComp+idxVs(2)) = vdq_s(2,:);
%
obj.x_v(rowIdxLoadComp+idxIabc(1)) = iabc(1,:);
obj.x_v(rowIdxLoadComp+idxIabc(2)) = iabc(2,:);
obj.x_v(rowIdxLoadComp+idxIabc(3)) = iabc(3,:);
%
obj.x_v(rowIdxLoadComp+idxWSlip) = wslip;
%
obj.x_v(rowIdxLoadComp+idxThetaSys) = thetaSys;
%
obj.x_v(rowIdxLoadComp+idxIs(1)) = idq_s(1,:);
obj.x_v(rowIdxLoadComp+idxIs(2)) = idq_s(2,:);
%
obj.x_v(rowIdxLoadComp+idxIr(1)) = idq_r(1,:);
obj.x_v(rowIdxLoadComp+idxIr(2)) = idq_r(2,:);
%
obj.x_v(rowIdxLoadComp+idxWr) = wr;
%
obj.x_v(rowIdxLoadComp+idxPQ(1)) = PQ_Zip(1,:);
obj.x_v(rowIdxLoadComp+idxPQ(2)) = PQ_Zip(2,:);
%
obj.x_v(rowIdxLoadComp+idxVmPreFilt) = Vm;
obj.x_v(rowIdxLoadComp+idxVm) = Vm;
%
obj.x_v(rowIdxLoadComp+idxIzip(1)) = ...
    real(loadComp_IComplexPhA_Zip);
obj.x_v(rowIdxLoadComp+idxIzip(2)) = ...
    imag(loadComp_IComplexPhA_Zip);

end