function [x1_0,x2_0,pref] = tgov1dValInitCalculations_assumeW0(obj,Pmech_0)
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

% Solve for pref
% (Laplace)
% (Pref - (Wr/Ws-1))*(1/R)*(1/ 1+s*T1) = X2
% (Pref - (Wr/Ws-1))*(1/R) = X2 + T1*s*X2
% (TD)
% pref - (wr/ws-1) = R*(x2 + T1*dx2/dt)
% (Steady-state derivatives 0, wr=ws)
% pref*(1/R) = x2
pref = obj.tgov1d_R.*x2_0;
end