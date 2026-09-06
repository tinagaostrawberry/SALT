function [Jvals,F,dParkV_dVphA,dFspeed_dWs] = stampOfIbrGFMMatrices_SS( ...
    ... % Current guess (pages: [nState x 1 x numIbrGFM])
    IphABC,vg_dq,vo_dq,if_dq,io_dq,vi_dq,gamma_dq,if_dq_ref,x_dq,wr, ...
    vo_dq_ref,pavg,qavg,theta,VphA, ...
    vo_dq_VR_avg,vo_dq_VR,vo_q_ref_droop, ws, ...
    ... % System constants (pages: [1 x 1 x numIbrGFM])
    Lf,rf,Cf,Lc,rc,kc_i,Gc,kv_i,Gv,Mp,Mq, pset,qset,vo_q_set, w_VR_rad,Rv, ...
    ... % System parameters
    w0, ibrGFMToBusI,ibrGFMToBusV)
%STAMPOFIBRGFMMATRICES_SS  Vectorised grid-forming IBR stamp, all units at once.
%
%   COMPUTATION:
%   J(x[k])*x[k+1] = J(x[k])*x[k]-F(x[k])
%
%   Jvals is [178 x numIbrGFM]; F is [33 x numIbrGFM].

n = size(theta,3);
Jvals = zeros(178,n);

% Pre-calc
stamp_j_invW0 = [0 -1/w0; 1/w0 0];      % stampReIm_2x2(1j/w0);
eye2  = eye(2);
eye2n = repmat(eye2,1,1,n);
stamp_j = [0 -1; 1 0];                  % stampReIm(1j);

% ----------------------------- Nonlinear ---------------------------------

% Park for I
% Io_phA = invT*Io_dq
cosTheta = cos(theta);
sinTheta = sin(theta);
invT  = [cosTheta, sinTheta; sinTheta,-cosTheta];   % [1,0;0,-1]*invT_gen
JinvT = [-sinTheta, cosTheta; cosTheta, sinTheta];
% Jacobian
dParkI_dIphA_dTheta_dIo = [-eye2n, ...
    pagemtimes(JinvT.*ibrGFMToBusI,io_dq), invT.*ibrGFMToBusI];
% RHS
FParkI = -IphABC(1:2,:,:) + pagemtimes(invT.*ibrGFMToBusI,io_dq);
Jvals(1:10,:) = reshape(dParkI_dIphA_dTheta_dIo,10,n);

% Ensure output current is balanced
% IphBC = phShiftBC*IphA
% phShiftBC = [stampReIm(exp(-1j*2*pi/3)); stampReIm(exp(1j*2*pi/3))];
phShiftBC = 0.5*[ ...
    -1, sqrt(3); -sqrt(3),-1; ...   % -120 deg
    -1,-sqrt(3);  sqrt(3),-1];      % +120 deg
% Jacobian
dFIphBC_dIphABC = [phShiftBC,-eye(4)];
% RHS
FIphBC = -IphABC(3:end,:,:) + pagemtimes(phShiftBC,IphABC(1:2,:,:));
Jvals(11:34,:) = repmat(dFIphBC_dIphABC(:),1,n);

% Park for V
% V_phA = invT*Vg_dq
% Jacobian
dParkV_dVphA = repmat(reshape(-eye2,[],1),1,n);
dParkV_dTheta_dVg = [pagemtimes(JinvT.*ibrGFMToBusV,vg_dq), invT.*ibrGFMToBusV];
% RHS
FParkV = -VphA + pagemtimes(invT.*ibrGFMToBusV,vg_dq);
Jvals(35:40,:) = reshape(dParkV_dTheta_dVg,6,n);

% Passive coupling
% io_dq' = (1/Lc)*(-rc*io_dq + vo_dq - vg_dq) + j*w*io_dq
% => ... = 0
jw = stamp_j_invW0.*wr;
% Jacobian
dFPassiveCoup_dIo_dWr_dVo_dVg = [(-rc./Lc).*eye2 + jw, ...
    pagemtimes(stamp_j_invW0,io_dq), (1./Lc).*eye2, -(1./Lc).*eye2];
% RHS
FPassiveCoup = (1./Lc).*(-rc.*io_dq + vo_dq - vg_dq) + pagemtimes(jw,io_dq);
Jvals(41:54,:) = reshape(dFPassiveCoup_dIo_dWr_dVo_dVg,14,n);

% Passive capacitor
% vo_dq' = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq + Rcap*(if_dq' - io_dq')
% =>
% 0 = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq
% Jacobian
dFPassiveCap_dIf_dIo_dVo_dW = [(1./Cf).*eye2, -(1./Cf).*eye2, jw, ...
    pagemtimes(stamp_j_invW0,vo_dq)];
% RHS
FPassiveCap = (1./Cf).*(if_dq - io_dq) + pagemtimes(jw,vo_dq);
Jvals(55:68,:) = reshape(dFPassiveCap_dIf_dIo_dVo_dW,14,n);

% Passive filter
% if_dq' = (1/Lf)*(-rf*if_dq + vi_dq - vo_dq) + j*w*if_dq
% => ... = 0
% Jacobian
dPassiveFilt_dIf_dWr_dVi_dVo = [(-rf./Lf).*eye2 + jw, ...
    pagemtimes(stamp_j_invW0,if_dq), (1./Lf).*eye2, -(1./Lf).*eye2];
% RHS
FPassiveFilt = (1./Lf).*(-rf.*if_dq + vi_dq - vo_dq) + pagemtimes(jw,if_dq);
Jvals(69:82,:) = reshape(dPassiveFilt_dIf_dWr_dVi_dVo,14,n);

% Current controller
% 1. vi_dq = kc_i*gamma_dq + kc_p*d(gamma_dq)/dt - j*wr*Lf*if_dq + Gc*vo_dq
% 2. d(gamma_dq)/dt = if_dq_ref - if_dq
% =>
% 1. vi_dq = kc_i*gamma_dq - j*wr*Lf*if_dq + Gc*vo_dq
% 2. 0 = if_dq_ref - if_dq
% Jacobian
dFCurrCont_dVi_dGamma_dIf_dWr_dVo = [-eye2n, kc_i.*eye2, ...
    -stamp_j_invW0.*(Lf.*wr), -Lf.*pagemtimes(stamp_j_invW0,if_dq), Gc.*eye2];
dFCurrIntErr_dIfRef_dIf = [eye2n,-eye2n];
% RHS
FibrCurrCont = -vi_dq + kc_i.*gamma_dq ...
    - Lf.*pagemtimes(jw,if_dq) + Gc.*vo_dq;
FibrCurrIntErr = if_dq_ref - if_dq;
Jvals(83:100,:)  = reshape(dFCurrCont_dVi_dGamma_dIf_dWr_dVo,18,n);
Jvals(101:108,:) = reshape(dFCurrIntErr_dIfRef_dIf,8,n);

% Voltage controller
% 1. if_dq_ref = kv_i*x_dq + kv_p*d(x_dq)/dt - j*wr*Cf*vo_dq + Gv*io_dq
% 2. d(x_dq)/dt = vo_dq_ref - vo_dq
% =>
% 1. if_dq_ref = kv_i*x_dq - j*wr*Cf*vo_dq + Gv*io_dq
% 2. 0 = vo_dq_ref - vo_dq
% Jacobian
dFVoltCont_dIfRef_dX_dVo_dWr_dIo = [-eye2n, kv_i.*eye2, ...
    -stamp_j_invW0.*(wr.*Cf), -Cf.*pagemtimes(stamp_j_invW0,vo_dq), Gv.*eye2];
dFVoltIntErr_dVoRef_dVo = [eye2n,-eye2n];
% RHS
FVoltCont = -if_dq_ref + kv_i.*x_dq - Cf.*pagemtimes(jw,vo_dq) + Gv.*io_dq;
FVoltIntErr = vo_dq_ref - vo_dq;
Jvals(109:126,:) = reshape(dFVoltCont_dIfRef_dX_dVo_dWr_dIo,18,n);
Jvals(127:134,:) = reshape(dFVoltIntErr_dVoRef_dVo,8,n);

% Droop
% 1. wr = w0 + Mp*(pset - pavg)
% 2. vo_q_ref_droop = vo_q_set + Mq*(qset-qavg)
% 3. pavg' = wmeas*(p-pavg) s.t. p = vg_dq.' * io_dq
%    NOTE: No factor of (3/2) is needed here. It would come from the
%    Park transform for peak, line-to-neutral, non-PU quantities; with
%    RMS, line-to-line, PU quantities it cancels.
% 4. qavg' = wmeas*(q-qavg) s.t. q = vg_dq.' * [0 -1; 1 0] * io_dq
% 5. d(theta)/dt = wr
% =>
% 1. wr = w0 + Mp*(pset - pavg)
% 2. vo_q_ref_droop = vo_q_set + Mq*(qset-qavg)
% 3. 0 = p - pavg
% 4. 0 = q - qavg
% 5. ws = wr
one = ones(1,1,n);
dFDroopW_dWr_dPavg          = -[one, Mp];
dFDroopV_dVoQRefDroop_dQavg = -[one, Mq];
dFDroopPavg_dPavg_dVg_dIo   = [-one, pagetranspose(io_dq), pagetranspose(vg_dq)];
dFDroopQavg_dQavg_dVg_Io    = [-one, ...
    pagemtimes(pagetranspose(io_dq),stamp_j.'), ...
    pagemtimes(pagetranspose(vg_dq),stamp_j)];
dFspeed_dWs = -ones(1,n);
% RHS
FDroopW    = -wr + w0 + Mp.*(pset - pavg);
FDroopV    = -vo_q_ref_droop + vo_q_set + Mq.*(qset - qavg);
FDroopPavg = pagemtimes(pagetranspose(vg_dq),io_dq) - pavg;
FDroopQavg = pagemtimes(pagemtimes(pagetranspose(vg_dq),stamp_j),io_dq) - qavg;
FSpeed     = -ws + wr;
Jvals(135:136,:) = reshape(dFDroopW_dWr_dPavg,2,n);
Jvals(137:138,:) = reshape(dFDroopV_dVoQRefDroop_dQavg,2,n);
Jvals(139:143,:) = reshape(dFDroopPavg_dPavg_dVg_dIo,5,n);
Jvals(144:148,:) = reshape(dFDroopQavg_dQavg_dVg_Io,5,n);
Jvals(149,:)     = 1;                       % dFspeed_dWr

% Virtual resistor
% ************************
% Vo_ref == vo_q_ref_droop
% ************************
% 1. d/dt(Vo_VR_avg) = w_VR*(Rv*Io - Vo_VR_avg)
% 2. Vo_VR = Vo_VR_avg - Rv*Io
% 3. Vo_ref = vo_ref_droop + Vo_VR
% =>
% 1. 0 = w_VR*(Rv*Io - Vo_VR_avg)
% 2. Vo_VR = Vo_VR_avg - Rv*Io
% 3. Vo_ref = vo_q_ref_droop + Vo_VR
w_VR = w_VR_rad/(2*pi);
dFVRVoAvg_dVoVrAvg_dIo    = w_VR.*[-eye2n, Rv.*eye2];
dFVRVo_dVoVR_dVoVrAvg_dIo = [-eye2n, eye2n, -Rv.*eye2];
dFVRVoRef_dVoRef_dVoVR    = [-eye2n, eye2n];
% RHS
FVRVoAvg  = w_VR.*(Rv.*io_dq - vo_dq_VR_avg);
FVRVo     = -vo_dq_VR + vo_dq_VR_avg - Rv.*io_dq;
FVRVoRef  = -vo_dq_ref + [zeros(1,1,n); vo_q_ref_droop] + vo_dq_VR;
Jvals(150:157,:) = reshape(dFVRVoAvg_dVoVrAvg_dIo,8,n);
Jvals(158:169,:) = reshape(dFVRVo_dVoVR_dVoVrAvg_dIo,12,n);
Jvals(170:177,:) = reshape(dFVRVoRef_dVoRef_dVoVR,8,n);
Jvals(178,:)     = 1;                       % dFVRVoRef_dVoQRefDroop

F = reshape([FParkI;FIphBC;FParkV; ...
    FPassiveCoup;FPassiveCap;FPassiveFilt; ...
    FibrCurrCont;FibrCurrIntErr;FVoltCont;FVoltIntErr; ...
    FDroopW;FDroopV;FDroopPavg;FDroopQavg; FSpeed; ...
    FVRVoAvg;FVRVo;FVRVoRef],33,n);
end
