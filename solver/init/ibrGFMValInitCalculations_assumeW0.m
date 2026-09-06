function  ...
    [vg_dq,vo_dq,vi_dq,if_dq,io_dq,gamma_dq,if_dq_ref,x_dq,wr, ...
    vo_dq_ref,pqavg,iabc_reIm, ...
    pset,qset,voDroop_q_set,THETA, ...
    vo_dq_VR_avg,vo_dq_VR,voDroop_q_ref] = ...
    ibrGFMValInitCalculations_assumeW0(obj, ...
    bus_VComplexPhA,bus_Va)

% Get bus
[~,busIdx_ibrGFM] = ismembertol(obj.ibrGFM_bus,obj.bus);

% Assumptions
w0 = obj.w0;
w0_PU = 1;

% Calc voltage and current (line-to-line (LL), RMS, bus PU):
ibrGFM_VComplexPhA_busPU = bus_VComplexPhA(busIdx_ibrGFM);
ibrGFM_SComplex_busPU = (obj.ibrGFM_MW + 1j*obj.ibrGFM_MVar)/ ...
    obj.baseMVA;
ibrGFM_IComplexPhA_busPU = ((ibrGFM_SComplex_busPU')./ ...
    (ibrGFM_VComplexPhA_busPU')).';
% Sanity check on current and power
assert(all(abs( ...
    abs(ibrGFM_IComplexPhA_busPU) ...
    - ...
    sqrt(obj.ibrGFM_MW.^2 + obj.ibrGFM_MVar.^2)/obj.baseMVA ...
    ./abs(ibrGFM_VComplexPhA_busPU) ...
    )<1e-6))

% Convert current to from bus to ibr PU
ibrGFM_IComplexPhA = ibrGFM_IComplexPhA_busPU.* ...
    obj.busToIbrGFM_I(busIdx_ibrGFM,1:obj.numIbrGFM);
ibrGFM_VComplexPhA = ibrGFM_VComplexPhA_busPU.* ...
    obj.busToIbrGFM_V(busIdx_ibrGFM,1:obj.numIbrGFM);

% Initialize DQ values
vg_dq = zeros(2,obj.numIbrGFM);
vo_dq = zeros(2,obj.numIbrGFM);
vi_dq = zeros(2,obj.numIbrGFM);
if_dq = zeros(2,obj.numIbrGFM);
io_dq = zeros(2,obj.numIbrGFM);
% Initialize controller values
gamma_dq = zeros(2,obj.numIbrGFM); % Current
if_dq_ref = zeros(2,obj.numIbrGFM);
x_dq = zeros(2,obj.numIbrGFM); % Voltage
% Initialize droop vars
wr = zeros(1,obj.numIbrGFM);
voDroop_q_ref = zeros(1,obj.numIbrGFM);
pqavg = zeros(2,obj.numIbrGFM);
% Initialize droop references (diff dim)
pset = zeros(obj.numIbrGFM,1);
qset = zeros(obj.numIbrGFM,1);
voDroop_q_set = zeros(obj.numIbrGFM,1);
% Initialize virtual resistor
vo_dq_VR_avg = zeros(2,obj.numIbrGFM);
vo_dq_VR = zeros(2,obj.numIbrGFM);
vo_dq_ref = zeros(2,obj.numIbrGFM);
% Initialize network
iabc_reIm = zeros(6,obj.numIbrGFM);
% Initialize theta
THETA = zeros(obj.numIbrGFM,1);
% Solve
for idxIbrGFM = 1:obj.numIbrGFM
    % NOTE: The underscore represents initial phase shift that enforces
    % d-axis voltage zero for vg; later will shift to enforce for vo as
    % convention

    % Park matrix
    THETA_ = bus_Va(busIdx_ibrGFM(idxIbrGFM)) - pi/2;
    T_ = getParkMatrices_2x3_PSCAD(THETA_);

    % Grid-interfacing
    [vgPhA,vgPhB,vgPhC] = calc3Ph(ibrGFM_VComplexPhA(idxIbrGFM));
    vg_dq_ = T_*real([vgPhA;vgPhB;vgPhC]);
    %
    [ioPhA,ioPhB,ioPhC] = calc3Ph(ibrGFM_IComplexPhA(idxIbrGFM));
    iabc_reIm(1:2:end,idxIbrGFM) = real([ioPhA;ioPhB;ioPhC])* ...
        obj.ibrGFMToBus_I(busIdx_ibrGFM(idxIbrGFM),idxIbrGFM);
    iabc_reIm(2:2:end,idxIbrGFM) = imag([ioPhA;ioPhB;ioPhC])* ...
        obj.ibrGFMToBus_I(busIdx_ibrGFM(idxIbrGFM),idxIbrGFM);
    io_dq_ = T_*real([ioPhA;ioPhB;ioPhC]);

    % Coupling
    % io_dq' = (1/Lc)*(-rc*io_dq + vo_dq - vg_dq) + j*w*io_dq
    % =>
    % vo_dq = vg_dq + rc*io_dq - j*w*Lc*io_dq
    vo_dq_ = vg_dq_ + obj.ibrGFM_Rc(idxIbrGFM)*io_dq_ ...
        - stampReIm(1j*w0_PU*obj.ibrGFM_Lc(idxIbrGFM))*io_dq_;

    % Perform phase shift to enforce vd is zero
    phShift = pi/2-angle(vo_dq_(1)+1j*vo_dq_(2));
    THETA(idxIbrGFM) = THETA_ + phShift;
    phShift_reIm = stampReIm(exp(1j*phShift));
    vo_dq(:,idxIbrGFM) = phShift_reIm*vo_dq_;
    vg_dq(:,idxIbrGFM) = phShift_reIm*vg_dq_;
    io_dq(:,idxIbrGFM) = phShift_reIm*io_dq_;
    % Sanity check that vo_dq is zero
    assert(abs(vo_dq(1,idxIbrGFM))<obj.absTolCheck)
    % Sanity check that coupling equation is still true
    assert(all(abs( ...
        vo_dq(:,idxIbrGFM) ...
        - ...
        (vg_dq(:,idxIbrGFM) + obj.ibrGFM_Rc(idxIbrGFM)*io_dq(:,idxIbrGFM) ...
        - stampReIm(1j*w0_PU*obj.ibrGFM_Lc(idxIbrGFM))*io_dq(:,idxIbrGFM)) ...
        )<obj.absTolCheck))
    % Sanity check that Park transform is still true
    T = getParkMatrices_2x3_PSCAD(THETA(idxIbrGFM));
    %invT = getinvParkMatrices_2x3_PSCAD(THETA(idxIbrGFM));
    assert(all(abs(vg_dq(:,idxIbrGFM)-T*real([vgPhA;vgPhB;vgPhC]))<obj.absTolCheck))
    assert(all(abs(io_dq(:,idxIbrGFM)-T*real([ioPhA;ioPhB;ioPhC]))<obj.absTolCheck))

    % Capacitor
    % vo_dq' = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq + Rcap*(if_dq' - io_dq')
    % =>
    % if_dq = io_dq - j*w*Cf*vo_dq
    if_dq(:,idxIbrGFM) = io_dq(:,idxIbrGFM) ...
        - stampReIm(1j*w0_PU*obj.ibrGFM_Cf(idxIbrGFM))*vo_dq(:,idxIbrGFM);

    % Filter
    % if_dq' = (1/Lf)*(-rf*if_dq + v_dq - vo_dq) + j*w*if_dq
    % =>
    % v_dq = vo_dq + rf*if_dq - j*w*Lf*if_dq
    vi_dq(:,idxIbrGFM) = vo_dq(:,idxIbrGFM) + obj.ibrGFM_Rf(idxIbrGFM)*if_dq(:,idxIbrGFM) ...
        - stampReIm(1j*w0_PU*obj.ibrGFM_Lf(idxIbrGFM))*if_dq(:,idxIbrGFM);

    % Current controller
    % 1. d(gamma_dq)/dt = if_dq_ref - if_dq
    % 2. vi_dq = kc_i*gamma_dq + kc_p*d(gamma_dq)/dt - j*wr*Lf*if_dq + Gc*vo_dq
    % =>
    % gamma_dq = (1/kc_i)*(vi_dq + j*wr*Lf*if_dq - Gc*vo_dq)
    gamma_dq(:,idxIbrGFM) = (1./obj.ibrGFM_kC_i(idxIbrGFM))*(vi_dq(:,idxIbrGFM) + ...
        stampReIm(1j*w0_PU*obj.ibrGFM_Lf(idxIbrGFM))*if_dq(:,idxIbrGFM) ...
        - obj.ibrGFM_GC(idxIbrGFM)*vo_dq(:,idxIbrGFM));
    % Sanity check
    if obj.ibrGFM_GC(idxIbrGFM)==1
        assert(all(abs(gamma_dq(:,idxIbrGFM)*obj.ibrGFM_kC_i(idxIbrGFM) ...
            -obj.ibrGFM_Rf(idxIbrGFM)*if_dq(:,idxIbrGFM))<obj.absTolCheck))
    end

    % Voltage controller
    % 1. d(x_dq)/dt = vo_dq_ref - vo_dq
    % 2. if_dq_ref = kv_i*x_dq + kv_p*d(x_dq)/dt - j*wr*Cf*vo_dq + Gv*io_dq
    % =>
    % x_dq = (1/kv_i)*(if_dq_ref + j*wr*Cf*vo_dq - Gv*io_dq)
    if_dq_ref(:,idxIbrGFM) = if_dq(:,idxIbrGFM);
    x_dq(:,idxIbrGFM) = (1/obj.ibrGFM_kV_i(idxIbrGFM))*(if_dq_ref(:,idxIbrGFM) + ...
        stampReIm(1j*w0_PU*obj.ibrGFM_Cf(idxIbrGFM))*vo_dq(:,idxIbrGFM) ...
        - obj.ibrGFM_GV(idxIbrGFM)*io_dq(:,idxIbrGFM));
    % Sanity check
    if obj.ibrGFM_GC(idxIbrGFM)==1
        assert(all(abs(x_dq(:,idxIbrGFM))<obj.absTolCheck))
    end

    % Droop
    % 1. wr = w + Mp*(pset - pavg)
    % 2. voDroop_q_ref = vo_q_set + Mq*(qset-qavg)
    % 3. pavg' = wmeas*(p-pavg) s.t. p = vo_dq.' * io_dq
    % 4. qavg' = wmeas*(q-qavg) s.t. q = [vd vq]*[0 -1; 1 0]*[id;iq]
    % 5. d(theta)/dt = wr
    wr(idxIbrGFM) = w0;
    voDroop_q_ref(idxIbrGFM) = vo_dq(2,idxIbrGFM);
    pqavg(:,idxIbrGFM) = [vg_dq(:,idxIbrGFM).'*io_dq(:,idxIbrGFM); ...
        vg_dq(:,idxIbrGFM).'*[0 -1; 1 0]*io_dq(:,idxIbrGFM)];
    % Control
    pset(idxIbrGFM) = pqavg(1,idxIbrGFM);
    qset(idxIbrGFM) = pqavg(2,idxIbrGFM);
    voDroop_q_set(idxIbrGFM) = vo_dq(2,idxIbrGFM);

    % Virtual resistor
    % d/dt(Vo_VR_avg) = w_VR*(Rv*Io - Vo_VR_avg)
    % =>
    % Vo_VR_avg = Rv*Io
    vo_dq_VR_avg(:,idxIbrGFM) = obj.ibrGFM_Rv(idxIbrGFM)*io_dq(:,idxIbrGFM);
    % Vo_VR = Vo_VR_avg - Rv*Io
    % =>
    % Vo_VR = 0
    vo_dq_VR(:,idxIbrGFM) = [0;0];
    % Vo_ref = vo_ref_droop + Vo_VR
    % =>
    % Vo_ref = [0;vo_ref_droop]
    vo_dq_ref(:,idxIbrGFM) = [0;voDroop_q_ref(idxIbrGFM)];
end

end
