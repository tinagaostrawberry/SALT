function [T,JT] = getParkMatrices_2x3_PSCAD(theta)
%GETPARKMATRICES_2X3_PSCAD  Power-invariant-free abc -> dq Park transform.
%
%   T  maps [a;b;c] to [d;q];  JT is dT/dtheta.

theta_phAbc = theta + [0 -2*pi/3 2*pi/3];
T  =  2/3*[ cos(theta_phAbc); sin(theta_phAbc)];
JT =  2/3*[-sin(theta_phAbc); cos(theta_phAbc)];

end
