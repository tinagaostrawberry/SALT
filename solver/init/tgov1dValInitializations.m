function [obj,pref_out] = tgov1dValInitializations(obj,tgov1dIndices)

% Check electrical and mechanical angles are equal
% theta = (p/2)*thetaM which means
% wr = (p/2)*wM, since w = dtheta/dt. If p = 2, then wr = (2/2)*wM = wM.
assert(obj.genrouP==2)
% When wr = wM, we have: Pmech := wM*TorqueM = wr*TorqueM = 1*TorqueM
%                                             |_|          |_|
%                                            (wr is normalized)
Pmech_0 = obj.genrou_TorqueM;

% ------------------------ Back propagate states --------------------------

% Backpropagate state x1
% (TD)
% x1 - (wr/ws-1)*Dt = pmech
% (Steady-state wr=ws)
% x1 = pmech
x1_0 = Pmech_0;

% Backpropagate state x2
% (Laplace)
% X1 = (1+s*T2)/(1+s*T3)*X2
% X1 + T3*s*X1 = X2 + T2*s*X2
% (TD)
% x1 + T3*dx1/dt = x2 + T2*dx2/dt
% (Steady-state derivatives 0)
% x2 = x1
x2_0 = x1_0;

% ------------------------- Back propagate pref ---------------------------

% Solve for pref
% (Laplace)
% (Pref - (Wr/Ws-1))*(1/R)*(1/ 1+s*T1) = X2
% (Pref - (Wr/Ws-1))*(1/R) = X2 + T1*s*X2
% (TD)
% pref - (wr/ws-1) = R*(x2 + T1*dx2/dt)
% (Steady-state derivatives 0, wr=ws)
% pref*(1/R) = x2
pref = obj.tgov1d_R.*x2_0;

% --------------------------------- Save ----------------------------------

% Default is saving for all governors
if nargin == 1
    tgov1dIndices = 1:obj.numTgov1d;
end

% Save states
idxRowTgov1d = obj.getRowIdx_fromTgov1dIndices(tgov1dIndices);
obj.x_v(idxRowTgov1d) = x1_0(tgov1dIndices);
obj.x_v(idxRowTgov1d+1) = x2_0(tgov1dIndices);
obj.x_v(idxRowTgov1d+2) = Pmech_0(tgov1dIndices);


% Output
pref_out = pref(tgov1dIndices);


end