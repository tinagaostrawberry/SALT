function errSmall = isErrSmall(x,x_prev,relTol,absTol)
%ISERRSMALL  Convergence test on a Newton step.
%
%   True when every entry satisfies either the relative or the absolute
%   tolerance, so entries near zero are not held to a relative test.

errSmall = all( (abs((x-x_prev)./x)<relTol) + (abs(x-x_prev)<absTol) ,'all');
end
