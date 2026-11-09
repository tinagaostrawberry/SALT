function [JvalsLin,F] = stampOfHbExcDc1Matrices_FD( ...
    ... % Current guess (pages: [nState x 1 x numExcDc1])
    x1,x2,x3,x4,x5,Vt, ...
    ... % Harmonic and frequency
    k,w0, ...
    ... % Parameters (pages: [1 x 1 x numExcDc1])
    Ke,Te,Tr,Tc,Tb,Ta,Ka,Tf1,Kf1, ...
    ... % Constants
    Vref)
%STAMPOFHBEXCDC1MATRICES_FD  Linear harmonic-balance EXDC1 stamp, one harmonic.
%
%   Harmonic-diagonal part of the EXDC1 equations, evaluated once per
%   harmonic k. Every input carries the exciters along the third dimension,
%   so the whole fleet is stamped in one pass (the legacy looped one device
%   at a time). JvalsLin is [44 x numExcDc1], one column per exciter, with
%   the eleven 2x2 (Re,Im) blocks flattened column-major and concatenated in
%   the fixed order documented below. F is the linear residual
%   [12 x numExcDc1] for this harmonic.
%
%   Ported from stampOfExcMatrices_linear_FD.m (legacy harmonic balance).
%   NOTE: the legacy also overwrote the imaginary DC rows of the device with
%   an identity to keep the matrix full rank. That is an assembly-level
%   overwrite, not a stamp, so it stays with the caller.

% VARIABLE LAYOUT (fixed): one EXDC1 block per harmonic is 12 rows, every
% quantity a (real,imag) pair. Within a device the row index is
% k*12 + localIndex, i.e. harmonic-outer / variable-inner - the ordering
% get_idxTermRow_ReIm.m builds in the legacy.
% INDICIES NOTATATION:
% Var name          Indices (Re,Im)
% x1                1:2
% x2                3:4
% x3                5:6
% x4                7:8
% x5                9:10
% Vt                11:12
idxHbExcDc1_x1 =  1:2;
idxHbExcDc1_x2 =  3:4;
idxHbExcDc1_x3 =  5:6;
idxHbExcDc1_x4 =  7:8;
idxHbExcDc1_x5 =  9:10;
idxHbExcDc1_Vt = 11:12;
nEqns = 12;

% Exciters run along dim 3. The linear part has no state dependence, so the
% current guess is only used to pin down the fleet size and to check that
% every page count agrees.
n = size(x1,3);
% Every parameter must be a page too, because the block values are reshaped
% straight out of the page dimension.
assert(all(cellfun(@(x)size(x,3)==n, ...
    {x2,x3,x4,x5,Vt,Ke,Te,Tr,Tc,Tb,Ta,Ka,Tf1,Kf1,Vref})))

% Harmonic frequency. k is an integer harmonic index, so the DC test is
% exact here and replaces the legacy's abs(wk)<1e-12 tolerance (the legacy
% received wk already multiplied out).
wk = k*w0;
isDc = (k==0);

% Block sizes, in the order they are concatenated into JvalsLin:
%  # | rows  | cols  | value
%  1 | 1:2   | 1:2   | stampReIm(-Ke - 1j*wk*Te)     ->  4
%  2 | 1:2   | 7:8   | stampReIm(1)                  ->  4
%  3 | 3:4   | 3:4   | stampReIm(1 + 1j*wk*Tr)       ->  4
%  4 | 3:4   | 11:12 | stampReIm(-1)                 ->  4
%  5 | 5:6   | 3:4   | stampReIm(1 + 1j*wk*Tc)       ->  4
%  6 | 5:6   | 5:6   | stampReIm(1 + 1j*wk*Tb)       ->  4
%  7 | 5:6   | 9:10  | stampReIm(1 + 1j*wk*Tc)       ->  4
%  8 | 7:8   | 5:6   | stampReIm(-Ka)                ->  4
%  9 | 7:8   | 7:8   | stampReIm(1 + 1j*wk*Ta)       ->  4
% 10 | 9:10  | 1:2   | stampReIm(-1j*wk*Kf1)         ->  4
% 11 | 9:10  | 9:10  | stampReIm(1 + 1j*wk*Tf1)      ->  4
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = [ ...
    numel(idxHbExcDc1_x1)^2; ...                              % 1
    numel(idxHbExcDc1_x1)*numel(idxHbExcDc1_x4); ...          % 2
    numel(idxHbExcDc1_x2)^2; ...                              % 3
    numel(idxHbExcDc1_x2)*numel(idxHbExcDc1_Vt); ...          % 4
    numel(idxHbExcDc1_x3)*numel(idxHbExcDc1_x2); ...          % 5
    numel(idxHbExcDc1_x3)^2; ...                              % 6
    numel(idxHbExcDc1_x3)*numel(idxHbExcDc1_x5); ...          % 7
    numel(idxHbExcDc1_x4)*numel(idxHbExcDc1_x3); ...          % 8
    numel(idxHbExcDc1_x4)^2; ...                              % 9
    numel(idxHbExcDc1_x5)*numel(idxHbExcDc1_x1); ...          % 10
    numel(idxHbExcDc1_x5)^2];                                 % 11
assert(isequal(nBlkVals.',repmat(4,1,11)))
% Running start offset of each block inside a JvalsLin column
offBlk = [0; cumsum(nBlkVals)];
nValsLin = offBlk(end);
assert(nValsLin==11*4)
assert(nValsLin==44)
JvalsLin = zeros(nValsLin,n);

% State 1 (backpropgation for x4)
% -Ke*x1(t) - Te*dx1(t)/dt + x4(t) = 0
% (Block 1) d/dx1
blk = stampReImPages(-Ke - 1j*wk*Te);
JvalsLin(offBlk(1)+(1:nBlkVals(1)),:) = reshape(blk,nBlkVals(1),n);
% (Block 2) d/dx4. Constant, so it is the same column for every exciter.
blk = stampReImPages(1);
JvalsLin(offBlk(2)+(1:nBlkVals(2)),:) = repmat(blk(:),1,n);

% State 2 (backpropagation for x2)
% x2(t) + Tr*dx2(t)/dt - et(t) = 0
% (Block 3) d/dx2
blk = stampReImPages(1 + 1j*wk*Tr);
JvalsLin(offBlk(3)+(1:nBlkVals(3)),:) = reshape(blk,nBlkVals(3),n);
% (Block 4) d/dVt. Constant.
blk = stampReImPages(-1);
JvalsLin(offBlk(4)+(1:nBlkVals(4)),:) = repmat(blk(:),1,n);

% State 3 (backpropagation for vref)
% x2(t) + Tc*dx2(t)/dt +
% x3(t) + Tb*dx3(t)/dt +
% x5(t) + Tc*dx5(t)/dt +
% - vref(t) - Tc*dvref(t)/dt = 0
%            |______________|
%                   0
% (Block 5) d/dx2
blk = stampReImPages(1 + 1j*wk*Tc);
JvalsLin(offBlk(5)+(1:nBlkVals(5)),:) = reshape(blk,nBlkVals(5),n);
% (Block 6) d/dx3
blk = stampReImPages(1 + 1j*wk*Tb);
JvalsLin(offBlk(6)+(1:nBlkVals(6)),:) = reshape(blk,nBlkVals(6),n);
% (Block 7) d/dx5. Same coefficient as block 5 (both carry Tc).
blk = stampReImPages(1 + 1j*wk*Tc);
JvalsLin(offBlk(7)+(1:nBlkVals(7)),:) = reshape(blk,nBlkVals(7),n);

% State 4 (backpropagation for x3)
% -Ka*x3(t) + x4(t) + Ta*dx4(t)/dt = 0
% (Block 8) d/dx3
blk = stampReImPages(-Ka);
JvalsLin(offBlk(8)+(1:nBlkVals(8)),:) = reshape(blk,nBlkVals(8),n);
% (Block 9) d/dx4
blk = stampReImPages(1 + 1j*wk*Ta);
JvalsLin(offBlk(9)+(1:nBlkVals(9)),:) = reshape(blk,nBlkVals(9),n);

% State 5 (backpropagation for x5)
% -Kf1*dx1(t)/dt + x5(t) + Tf1*dx5(t)/dt = 0
% (Block 10) d/dx1. Exactly zero at DC, so one branchless expression covers
% both the DC and the higher-harmonic case.
blk = stampReImPages(-1j*wk*Kf1);
JvalsLin(offBlk(10)+(1:nBlkVals(10)),:) = reshape(blk,nBlkVals(10),n);
% (Block 11) d/dx5
blk = stampReImPages(1 + 1j*wk*Tf1);
JvalsLin(offBlk(11)+(1:nBlkVals(11)),:) = reshape(blk,nBlkVals(11),n);

% ------------------------- Linear residual -------------------------------
% RESIDUAL SIGN. The legacy assembles Y*x = I. The caller here works with
% J*x - F, i.e. J*x = F, and J is the legacy's Yexc, so F is the legacy's
% Iexc with NO sign change. The only linear term with no state attached is
% the "- vref(t)" of state 3, which moves to the right-hand side as +Vref -
% exactly the legacy's "Iexc(5) = Vref" under "if isDc".
% It lands on the real row of state 3 (row 5) and only at k=0 because vref
% is a pure DC quantity. No N scaling appears here: the caller already
% multiplies vref by N when it builds exc_vref0, so the value arriving in
% Vref is in the same scaled-DFT-coefficient convention as the state vector.
F = zeros(nEqns,n);
if isDc
    F(5,:) = reshape(Vref,1,n); % Only real
end

end


function Y_ReIm = stampReImPages(Y)
%STAMPREIMPAGES  Page-aware stampReIm, i.e. [real(Y) -imag(Y); imag(Y) real(Y)]
%   written with the real and imaginary rows/cols interleaved. Identical to
%   stampReIm.m, extended with a trailing page dimension so one call covers
%   every exciter.

[nr,nc,n] = size(Y);
Y_ReIm = zeros(2*nr,2*nc,n);

Y_ReIm(1:2:end,1:2:end,:) =  real(Y);
Y_ReIm(1:2:end,2:2:end,:) = -imag(Y);
Y_ReIm(2:2:end,1:2:end,:) =  imag(Y);
Y_ReIm(2:2:end,2:2:end,:) =  real(Y);

end
