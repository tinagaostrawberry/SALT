function [Jvals,Jwr,Ilin, FTorqueM,JTorqueM_dPmech_dWr_dTm] = ...
    stampOfTgov1dMatrices_SS( ...
    ... % Current guess (rows are devices)
    pmech,Tm,wr, ...
    ... % System constants
    R,Dt,pref, ...
    ... % System parameters
    w0, genrouP)
%STAMPOFTGOV1DMATRICES_SS  Vectorised TGOV1D stamp for every governor at once.
%
%   Jvals [5 x numTgov1d]  in-block entries
%   Jwr   [3 x numTgov1d]  coupling of each governor row to machine speed
%   Ilin  [3 x numTgov1d]  constant right-hand side
%   The last two outputs drive the machine's mechanical-torque equation.

n = numel(pmech);
Jvals = zeros(5,n);
Ilin  = zeros(3,n);

% ------------------------------ Linear -----------------------------------

% State 1
% x1 + T3*dx1/dt - x2 - T2*dx2/dt = 0
% =>
% x1 - x2 = 0
% Linear stamp (w.r.t x1,x2)
Jvals(1,:) =  1;
Jvals(2,:) = -1;

% State 2
% pref - (wr/w0-1) - R*x2 - R*T1*dx2/dt = 0
% =>
% pref - (wr/w0-1) - R*x2 = 0
% Linear stamp (w.r.t x2,wr)
Jvals(3,:) = -R(:).';
% Linear RHS
Ilin(2,:) = -pref(:).' - 1;

% State 3
% x1 - (wr/w0-1)*Dt - pmech = 0
% Linear stamp (w.r.t x1,x3/pmech,wr)
Jvals(4,:) =  1;
Jvals(5,:) = -1;
% Linear RHS
Ilin(3,:) = -Dt(:).';

% Coupling to machine speed, one entry per governor row
Jwr = [zeros(1,n); -ones(1,n)/w0; -Dt(:).'/w0];

% ----------------------------- Nonlinear ---------------------------------
% COMPUTATION:
% J(x[k])*x[k+1] = J(x[k])*x[k]-F(x[k])

% Tm
% Check electrical and mechanical angles are equal
% theta = (p/2)*thetaM which means
% wr = (p/2)*wM, since w = dtheta/dt. If p = 2, then wr = (2/2)*wM = wM.
assert(genrouP==2)
% When wr = wM, we have:
% pmech - (wr/w0)*Tm = 0
% Jacobian
JTorqueM_dPmech_dWr_dTm = [ones(1,n); -Tm(:).'/w0; -wr(:).'/w0];
% RHS
FTorqueM = pmech(:).' - (wr(:).'/w0).*Tm(:).';
end
