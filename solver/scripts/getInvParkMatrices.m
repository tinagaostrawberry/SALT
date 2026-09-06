function [invT,JinvT] = getInvParkMatrices(theta)
%GETINVPARKMATRICES  dq0 -> abc inverse Park transform.
%
%   invT maps [d;q;0] to [a;b;c];  JinvT is d(invT)/dtheta, with
%   d(theta)/dt = wr.

theta_phAbc = theta + [0;-2*pi/3;2*pi/3];
invT  = [ cos(theta_phAbc) -sin(theta_phAbc) ones(3,1)];
JinvT = [-sin(theta_phAbc) -cos(theta_phAbc) zeros(3,1)];

end
