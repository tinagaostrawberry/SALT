function [x1,x2,x3,x4,x5,vref] = hbExcDc1ValInitCalculations(obj,et0,efd)
%HBEXCDC1VALINITCALCULATIONS  Back-propagate EXDC1 states at steady state.
%
%   et0 is the machine terminal voltage, efd its field voltage, both as
%   TIME-DOMAIN values. Every derivative is zero at steady state, so each
%   block inverts to an algebraic relation and the chain can be walked
%   backwards from efd to vref.

exc_efd = efd.*obj.excDc1_genToExc;

% x1: efd = x1*wr_pu, and wr_pu = 1 at steady state
x1 = exc_efd;

% x4: Te*dx1/dt = x4 - Ke*x1  =>  x4 = Ke*x1
x4 = obj.excDc1_Ke.*x1;

% x3: x4 + Ta*dx4/dt = Ka*x3  =>  x3 = x4/Ka
x3 = x4./obj.excDc1_Ka;

% x5: x5 + Tf1*dx5/dt = Kf1*dx1/dt  =>  x5 = 0
x5 = zeros(obj.numExcDc1,1);

% x2: x2 + Tr*dx2/dt = et  =>  x2 = et
x2 = et0;

% vref: x3 + Tb*dx3/dt = (vref-x2-x5) + Tc*d(vref-x2-x5)/dt
%   =>  vref = x3 + x2 + x5
vref = x3 + x2 + x5;

end
