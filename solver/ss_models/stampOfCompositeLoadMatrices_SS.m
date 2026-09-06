function [Jvals,F, dFParkV_dVphA,dFslip_dwSys,dFs_dwSys,dFswing_dwSys] = ...
    stampOfCompositeLoadMatrices_SS( ...
    ... % States (pages: [nState x 1 x numLoadComp])
    vs,iphABC,wSlip,thetaSys,is,ir,wr,vphA,wSys, PQZip,Vm,iZip, ...
    ... % Parameters
    R,L,L_Te,T0,m, PQz,PQi,PQp, ...
    ... % Consts
    w0, busToLoadCompI,busToLoadCompV)
%STAMPOFCOMPOSITELOADMATRICES_SS  Vectorised IM + ZIP composite-load stamp.
%
%   Jvals is [99 x numLoadComp]; F is [20 x numLoadComp]. The trailing
%   outputs are the couplings into the terminal bus voltage and into the
%   system-frequency column.

n = size(wr,3);
Jvals = zeros(99,n);
eye2  = eye(2);
eye2n = repmat(eye2,1,1,n);
one   = ones(1,1,n);

% ------------------------- Induction motor -------------------------------

% Park V (eqns 1-2)
% vdq_s = T(theta_sys)*vabc
% =>
% vdq_s = vphA
% Jacobian
dFParkV_dVs   = -eye2n;
dFParkV_dVphA = reshape(eye2.*busToLoadCompV,4,n);
% RHS
FParkV = -vs + vphA.*busToLoadCompV;
Jvals(1:4,:) = reshape(dFParkV_dVs,4,n);

% Park I (eqns 3-4)
% (idq_s+idq_ZIP) = T(theta_sys)*iabc
% =>
% (idq_s+idq_ZIP) = iphA
% Jacobian
dFParkI_dIs_dIzip_dIphA = [-eye2n,-eye2n, eye2.*busToLoadCompI];
% RHS
FParkI = -(is + iZip) + iphABC(1:2,:,:).*busToLoadCompI;
Jvals(5:16,:) = reshape(dFParkI_dIs_dIzip_dIphA,12,n);

% Ensure output current is balanced (eqns 5-8)
% IphBC = phShiftBC*IphA
% s.t. phShiftBC = [stampReIm(exp(-1j*2*pi/3)); stampReIm(exp(1j*2*pi/3))]
phShiftBC = 0.5*[ ...
    -1, sqrt(3); -sqrt(3),-1; ...   % -120 deg
    -1,-sqrt(3);  sqrt(3),-1];      % +120 deg
% Jacobian
dFIphBC_dIphABC = [phShiftBC,-eye(4)];
% RHS
FIphBC = -iphABC(3:end,:,:) + pagemtimes(phShiftBC,iphABC(1:2,:,:));
Jvals(17:40,:) = repmat(dFIphBC_dIphABC(:),1,n);

% Slip frequency (eqn 9)
% wSlip = wSys - wr
% Jacobian
dFslip_dwSys = ones(1,n);
% RHS
Fslip = -wSlip + wSys - wr;
Jvals(41:42,:) = -1;                    % dFslip_dwSlip_dwr

% System angle (eqn 10)
% d(thetasys)/dt = wSys  =>  thetasys = wSys*t  =>  thetasys = 0 @ t=0
% Jacobian
% RHS
Fsys = -thetaSys;
Jvals(43,:) = -1;                       % dFsys_dThetaSys

% Stator/rotor equations (eqn 11-14)
% |vs| = [R]*|idq_s| + (d/dt)*[L]/w0*|idq_s| + [wMat]/w0*[L]*|idq_s|
% |0 |       |idq_r|                 |idq_r|                 |idq_r|
% =>
% |vs| = [R]*|idq_s| + [wMat]/w0*[L]*|idq_s|
% |0 |       |idq_r|                 |idq_r|
%
% s.t.
% R = |Rs_I2x2      |
%     |      Rr_I2x2|
%
% L = |Lss_I2x2 Lm_I2x2 |
%     |Lm_I2x2  Lrr_I2x2|
vsr = [vs; zeros(2,1,n)];
isr = [is; ir];
wMat = zeros(4,4,n);
wMat(1,2,:) = -wSys;  wMat(2,1,:) =  wSys;
wMat(3,4,:) = -wSlip; wMat(4,3,:) =  wSlip;
L_PU  = L/w0;
phisr = pagemtimes(L_PU,isr);
% Jacobian
dFsr_dvs_disr_dwSlip = [ ...
    repmat([-eye2;zeros(2)],1,1,n), ...
    R + pagemtimes(wMat,L_PU), ...
    [zeros(2,1,n); -phisr(4,:,:); phisr(3,:,:)]];
dFs_dwSys = reshape([-phisr(2,:,:); phisr(1,:,:)],2,n);
% RHS
Fsr = -vsr + pagemtimes(R,isr) + pagemtimes(wMat,phisr);
Jvals(44:71,:) = reshape(dFsr_dvs_disr_dwSlip,28,n);

% Swing equation (eqn 15)
% (2*H/w0)*d(wr)/dt = ((isr^T)*L_Te*isr) - T0*(wr/wsys)^m
% =>
% 0 = ((isr^T)*L_Te*isr) - T0*(wr/wsys)^m
LTeSym = L_Te + pagetranspose(L_Te);
% Jacobian
dFswing_dwr_disr = [-T0.*(m.*wr.^(m-1))./(wSys.^m), ...
    pagemtimes(pagetranspose(isr),LTeSym)];
dFswing_dwSys = reshape(-T0.*wr.^m.*(-m).*wSys.^(-m-1),1,n);
% RHS
Fswing = pagemtimes(pagemtimes(pagetranspose(isr),L_Te),isr) ...
    - T0.*(wr./wSys).^m;
Jvals(72:76,:) = reshape(dFswing_dwr_disr,5,n);

% ---------------------------- ZIP model ----------------------------------

% ZIP (eqns 16-17)
% PQ = PQz*Vm^2 + PQi*Vm + PQp
% Jacobian
dFZip_dPQ_dVm = [-eye2n, 2*PQz.*Vm + PQi];
% RHS
FZip = -PQZip + PQz.*Vm.^2 + PQi.*Vm + PQp;
Jvals(77:82,:) = reshape(dFZip_dPQ_dVm,6,n);

% Vm (eqn 18)
% Vm = (vd_s^2+vq_s^2)^(1/2)
vsMag = sqrt(vs(1,:,:).^2 + vs(2,:,:).^2);
% Jacobian
dFVm_dvm_dvs = [-one, pagetranspose(vs)./vsMag];
% RHS
FVm = -Vm + vsMag;
Jvals(83:85,:) = reshape(dFVm_dvm_dvs,3,n);

% Izip (eqn 19-20)
% idZip = (P*vd + Q*vq)/Vm^2
% iqZip = (P*vq - Q*vd)/Vm^2
iZipEval = pagemtimes([vs, [vs(2,:,:); -vs(1,:,:)]],PQZip);
Vm2_inv  = 1./Vm.^2;
% Jacobian
blk = [pagetranspose(vs), pagetranspose(PQZip); ...
       vs(2,:,:), -vs(1,:,:), -PQZip(2,:,:), PQZip(1,:,:)];
dFIZip_dPQ_dvs_dVm_diZip = ...
    [[blk, -2*iZipEval./Vm].*Vm2_inv, -eye2n];
% RHS
FIZip = -iZip + iZipEval.*Vm2_inv;
Jvals(86:99,:) = reshape(dFIZip_dPQ_dvs_dVm_diZip,14,n);

F = reshape([FParkV;FParkI;FIphBC;Fslip;Fsys;Fsr;Fswing;FZip;FVm;FIZip],20,n);
end
