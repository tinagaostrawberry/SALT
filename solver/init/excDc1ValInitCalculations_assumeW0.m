function [x1_0,x2_0,x3_0,x4_0,x5_0,vref] = ...
    excDc1ValInitCalculations_assumeW0(obj,et_0,genrou_efd)
excDc1_efd = genrou_efd.*obj.excDc1_genrouToExcDc1;

% Back propagate state x1
% x1 = efd ./ w_pu = efd / 1 = efd
x1_0 = excDc1_efd;

% Back propagate state x4
% (Laplace)
% X1(s) = ( 1 / s*Te )*(X4(s)-Ke*X1(s))
% (TD)
% Te*dx1(t)/dt = x4(t) - Ke*x1(t)
% (Steady-state derivatives 0)
% x4(t) = Ke * x1(t)
x4_0 = obj.excDc1_Ke.*x1_0;

% Back propagate state x3
% (Laplace)
% X4(s) = ( Ka / 1+s*Ta )*X3(s)
% (TD)
% x4(t) + Ta*dx4(t)/dt = Ka*x3(t)
% (Steady-state derivatives 0)
% x3(t) = x4(t)/Ka
x3_0 = x4_0./obj.excDc1_Ka;

% Back propagate state x5
% (Laplace)
% X5(s) = ( s*Kf1 / 1+s*Tf1 )*X1(s)
% (TD)
% x5(t) + Tf1*dx5(t)/dt = Kf1*dx1(t)/dt
% (Steady-state derivatives 0)
% x5(t) = 0
x5_0 = zeros(obj.numExcDc1,1);

% Back propagate state x2
% (Laplace)
% X2(s) = (1 / 1+s*Tr )*Et(s)
% (TD)
% x2(t) + Tr*dx2(t)/dt = et(t)
% (Steady-state derivatives 0)
% x2(t) = et(t)
x2_0 = et_0;

% Solve for vref
% (Laplace)
% X3(s) = ( 1+s*Tc / 1+s*Tb )*(Vref(s)-X2(s)-X5(s))
% (TD)
% x3(t) + Tb*dx3(t)/dt = (vref(t)-x2(t)-x5(t)) + Tc*d(vref(t)-x2(t)-x5(t))/dt
% (Steady-state derivatives 0)
% vref(t) = x3(t) + x2(t) + x5(t)
vref = x3_0 + x2_0 + x5_0;
end