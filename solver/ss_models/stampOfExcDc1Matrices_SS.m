function [Jvals,Ilin, FexcDc1Vt,JexcDc1Vt_dvt,JexcDc1Vt_dvabc, ...
    FexcDc1Efd,JexcDc1Efd_x1,JexcDc1Efd_efd,JexcDc1Efd_wr] = ...
    stampOfExcDc1Matrices_SS( ...
    ... % Current guess (rows are devices)
    x1,genrou_efd,vt_genrouPU,vtabc_busPU,wr, ...
    ... % System constants
    Ka,Ke,vref, ...
    ... % System parameters
    w0, busToGenrouV,excDc1_genrouToExcDc1)
%STAMPOFEXCDC1MATRICES_SS  Vectorised EXDC1 stamp for every exciter at once.
%
%   Jvals is [10 x numExcDc1], holding the in-block entries in the order of
%   the triplet map. The remaining outputs are the coupling terms into the
%   terminal-voltage row and into the machine's field-voltage equation.
%   vtabc_busPU is [3 x numExcDc1].

n = numel(x1);
Jvals = zeros(10,n);
Ilin  = zeros(5,n);

% ------------------------------- Linear ----------------------------------

% State 1
% -Ke*x1(t) - Te*dx1(t)/dt + x4(t) = 0
% =>
% -Ke*x1(t) + x4(t) = 0
% Linear Ystamp (w.r.t x1,x4)
Jvals(1,:) = -Ke(:).';
Jvals(2,:) = 1;

% State 2
% x2(t) + Tr*dx2(t)/dt - vt(t) = 0
% =>
% x2(t) - vt(t) = 0
% Linear Ystamp (w.r.t x2,vt)
Jvals(3,:) =  1;
Jvals(4,:) = -1;

% State 3
% x2(t) + Tc*dx2(t)/dt +
% x3(t) + Tb*dx3(t)/dt +
% x5(t) + Tc*dx5(t)/dt +
% - vref(t) - Tc*dvref(t)/dt = 0
%            |______________|
%                   0
% =>
% x2(t) + x3(t) + x5(t) - vref(t) = 0
% Linear Ystamp (w.r.t x2,x3,x5)
Jvals(5:7,:) = 1;
% Linear RHS
Ilin(3,:) = vref(:).';

% State 4
% -Ka*x3(t) + x4(t) + Ta*dx4(t)/dt = 0
% =>
% -Ka*x3(t) + x4(t) = 0
% Linear Ystamp (w.r.t x3,x4)
Jvals(8,:) = -Ka(:).';
Jvals(9,:) = 1;

% State 5
% -Kf1*dx1(t)/dt + x5(t) + Tf1*dx5(t)/dt = 0
% =>
% x5(t) = 0
% Linear Ystamp (w.r.t x5)
Jvals(10,:) = 1;

% ----------------------------- Nonlinear ---------------------------------
% COMPUTATION:
% J(x[k])*x[k+1] = J(x[k])*x[k]-F(x[k])

% Efd normalized by speed
% exc_efd - x1*wr_pu = 0
% exc_efd = genToExc*gen_efd
% =>
% genToExc*gen_efd - x1*wr_pu = 0
% Jacobian
JexcDc1Efd_x1  = -wr(:).'/w0;
JexcDc1Efd_efd = excDc1_genrouToExcDc1(:).';
JexcDc1Efd_wr  = -x1(:).'/w0;
% RHS
FexcDc1Efd = (genrou_efd(:).').*JexcDc1Efd_efd - (x1(:).').*wr(:).'/w0;

% Vt
% -Vt + sqrt(2/3)*(va^2 + vb^2 + vc^2)^(1/2) = 0 s.t. va,vb,vc are gen PU
% Jacobian
c = busToGenrouV(:).';                  % 1 x n
vtabc_genrouPU = c.*vtabc_busPU;        % 3 x n (real part only)
sumSq = sum(vtabc_genrouPU.^2,1);       % 1 x n
JexcDc1Vt_dvt = -ones(1,n);
% The trailing conversion is for when the derivative multiplies vabc
JexcDc1Vt_dvabc = sqrt(2/3)*(sumSq.^(-1/2)).*vtabc_genrouPU.*c;   % 3 x n
% RHS
FexcDc1Vt = -vt_genrouPU(:).' + sqrt(2/3)*(sumSq.^(1/2));
end
