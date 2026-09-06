function [Jvals,F,dParkV_dVphA,dFspeed_dWs] = stampOfGenrouMatrices_SS( ...
    ... % Current guess (pages: [nState x 1 x numGenrou])
    i,e,IphABC,VphA,theta,wr,Tm, ws, ...
    ... % System constants
    Rmat,Lspmat, ...
    ... % System parameters
    w0, genrouToBusI,genrouToBusV)
%STAMPOFGENROUMATRICES_SS  Vectorised GENROU stamp for every machine at once.
%
%   Every input carries the machines along the third dimension, so the
%   whole fleet is stamped in one pass. Jvals is [68 x numGenrou], one
%   column per machine, ordered to match the triplet map built by
%   SALT.buildStampMaps. F is the residual [15 x numGenrou].

n = size(theta,3);
Jvals = zeros(68,n);

% Stator                                [eqns 1-3]
% E = -[R]*I + [L]*j*w*I + [Lsp]*I*(wr/w0)
%              |_______|
%                  0
idxDqf = [1 2 4];
R   = Rmat(idxDqf,idxDqf,:);
Lsp = Lspmat(idxDqf,idxDqf,:);
A = -R + Lsp.*(wr/w0);
% Jacobian
dFstator_dE_dI_dWr = [repmat(-eye(3),1,1,n), A, pagemtimes(Lsp,i)/w0];
% RHS
Fstator = -e + pagemtimes(A,i);
% Save Jacobian values
Jvals(1:21,:) = reshape(dFstator_dE_dI_dWr,21,n);

% Park transform for voltage                  [eqns 4-5]
% VphA = invT*Edq
cosTheta = cos(theta);
sinTheta = sin(theta);
invT  = [cosTheta,-sinTheta; sinTheta, cosTheta];
JinvT = [-sinTheta,-cosTheta; cosTheta,-sinTheta];
% Jacobian
dParkV_dVphA = repmat(reshape(-eye(2),[],1),1,n);
dParkV_dTheta_dEdq = [pagemtimes(JinvT.*genrouToBusV,e(1:2,:,:)), ...
    invT.*genrouToBusV];
% RHS
FParkV = -VphA + pagemtimes(invT.*genrouToBusV,e(1:2,:,:));
% Save Jacobian values
Jvals(22:27,:) = reshape(dParkV_dTheta_dEdq,6,n);

% Efd is from the exciter, so leave empty! Stamped in the exciter equations

% Rotor speed is equal to system speed  [eqn 7]
% wr = ws
% Jacobian
Jvals(28,:) = -1;                       % dFspeed_dWr
dFspeed_dWs = ones(1,n);
% RHS
Fspeed = -wr + ws;

% Swing  [eqn 8]
% Tm = Pt + Ra*It^2 = (ed*id+eq*iq) + Ra*(id^2+iq^2)
% ...but when wr=/=w0, need normalization because derivation of electric
% torque uses speed terms from the stator equation
% =>
% Tm*(wr/w0) = (ed*id+eq*iq) + Ra*(id^2+iq^2)
Ra  = Rmat(1,1,:);
idq = i(1:2,:,:);
edq = e(1:2,:,:);
% Jacobian
dFswing_dTm_dWr_dEdq_dIdq = [-wr/w0, -Tm/w0, pagetranspose(idq), ...
    pagetranspose(edq) + 2*Ra.*pagetranspose(idq)];
% RHS
Fswing = -Tm.*(wr/w0) + pagemtimes(pagetranspose(edq),idq) ...
    + Ra.*pagemtimes(pagetranspose(idq),idq);
% Save Jacobian values
Jvals(29:34,:) = reshape(dFswing_dTm_dWr_dEdq_dIdq,6,n);

% Park transform for current  [eqns 9-10]
% IphA = invT*Idq
% Jacobian
dParkI_dIphA_dTheta_dIdq = [repmat(-eye(2),1,1,n), ...
    pagemtimes(JinvT.*genrouToBusI,idq), invT.*genrouToBusI];
% RHS
FParkI = -IphABC(1:2,:,:) + pagemtimes(invT.*genrouToBusI,idq);
% Save Jacobian values
Jvals(35:44,:) = reshape(dParkI_dIphA_dTheta_dIdq,10,n);

% Ensure output current is balanced  [eqns 11-14]
% IphBC = phShiftBC*IphA
% phShiftBC = [stampReIm(exp(-1j*2*pi/3)); stampReIm(exp(1j*2*pi/3))];
phShiftBC = 0.5*[ ...
    -1, sqrt(3); -sqrt(3),-1; ...   % -120 deg
    -1,-sqrt(3);  sqrt(3),-1];      % +120 deg
% Jacobian (constant)
dFIphBC_dIphABC = [phShiftBC,-eye(4)];
% RHS
FIphBC = -IphABC(3:end,:,:) + pagemtimes(phShiftBC,IphABC(1:2,:,:));
% Save Jacobian values
Jvals(45:68,:) = repmat(dFIphBC_dIphABC(:),1,n);

% Tm is from the governor, so leave empty! Stamped in the gov equations

F = reshape([Fstator;FParkV;zeros(1,1,n);Fspeed;Fswing;FParkI;FIphBC; ...
    zeros(1,1,n)],15,n);
end
