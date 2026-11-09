function [x1,x2,pmech,pref] = hbTgov1dValInitCalculations(obj,TorqueM)
%HBTGOV1DVALINITCALCULATIONS  Back-propagate TGOV1D states at steady state.
%
%   TorqueM is the machine's mechanical torque as a TIME-DOMAIN value.
%
%   theta = (p/2)*thetaM (Kundur 3.1.2), so wr = (p/2)*wM. With p = 2 that
%   is wr = wM, which makes pmech := wM*TorqueM = wr*TorqueM = TorqueM once
%   wr is normalised. The assert is what keeps that substitution honest.
assert(obj.genP==2)

pmech = TorqueM;

% x1: x1 - (wr/ws-1)*Dt = pmech, and wr = ws at steady state
x1 = pmech;

% x2: x1 + T3*dx1/dt = x2 + T2*dx2/dt  =>  x2 = x1
x2 = x1;

% pref: (pref - (wr/ws-1))/R = x2 + T1*dx2/dt  =>  pref = R*x2
pref = obj.tgov1d_R.*x2;

end
