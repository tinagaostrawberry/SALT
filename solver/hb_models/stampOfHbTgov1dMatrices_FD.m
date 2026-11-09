function [JvalsLin,F,dFx2_dWr,dFpmech_dWr] = stampOfHbTgov1dMatrices_FD( ...
    ... % Current guess (pages: [nState x 1 x numTgov1d])
    x1,x2,Pmech, ...
    ... % Harmonic and frequency
    k,w0, ...
    ... % Parameters (pages: [1 x 1 x numTgov1d])
    T1,T2,T3,R,Dt, ...
    ... % Constants
    w0nom,N,Pref)
%STAMPOFHBTGOV1DMATRICES_FD  Linear harmonic-balance TGOV1D stamp, one harmonic.
%
%   Harmonic-diagonal part of the TGOV1D equations, evaluated once per
%   harmonic k. Every input carries the governors along the third dimension,
%   so the whole fleet is stamped in one pass (the legacy looped one device
%   at a time). JvalsLin is [20 x numTgov1d], one column per governor, with
%   the five 2x2 (Re,Im) blocks flattened column-major and concatenated in
%   the fixed order documented below. F is the linear residual
%   [6 x numTgov1d] for this harmonic. dFx2_dWr and dFpmech_dWr are each
%   [4 x numTgov1d]: the couplings out of the device block into the GENROU
%   speed Wr at the SAME harmonic (the legacy's Ygov_wr).
%
%   Ported from stampOfGovMatrices_linear_FD.m (legacy harmonic balance).
%   NOTE: the legacy also overwrote the imaginary DC rows of the device with
%   an identity to keep the matrix full rank. That is an assembly-level
%   overwrite, not a stamp, so it stays with the caller.

% VARIABLE LAYOUT (fixed): one TGOV1D block per harmonic is 6 rows, every
% quantity a (real,imag) pair. Within a device the row index is
% k*6 + localIndex, i.e. harmonic-outer / variable-inner - the ordering
% get_idxTermRow_ReIm.m builds in the legacy.
% INDICIES NOTATATION:
% Var name          Indices (Re,Im)
% x1                1:2
% x2                3:4
% Pmech             5:6
idxHbTgov1d_x1    = 1:2;
idxHbTgov1d_x2    = 3:4;
idxHbTgov1d_Pmech = 5:6;
nEqns = 6;

% Governors run along dim 3. The linear part has no state dependence, so the
% current guess is only used to pin down the fleet size and to check that
% every page count agrees.
n = size(x1,3);
% Every parameter must be a page too, because the block values are reshaped
% straight out of the page dimension.
assert(all(cellfun(@(x)size(x,3)==n, ...
    {x2,Pmech,T1,T2,T3,R,Dt,Pref})))

% Harmonic frequency. k is an integer harmonic index, so the DC test is
% exact here and replaces the legacy's abs(wk)<1e-12 tolerance (the legacy
% received wk already multiplied out).
wk = k*w0;
isDc = (k==0);

% Block sizes, in the order they are concatenated into JvalsLin:
%  # | rows | cols | value
%  1 | 1:2  | 1:2  | stampReIm(1 + 1j*wk*T3)      -> 4
%  2 | 1:2  | 3:4  | stampReIm(-1 - 1j*wk*T2)     -> 4
%  3 | 3:4  | 3:4  | stampReIm(-R - 1j*wk*R*T1)   -> 4
%  4 | 5:6  | 1:2  | stampReIm(1)                 -> 4
%  5 | 5:6  | 5:6  | stampReIm(-1)                -> 4
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = [ ...
    numel(idxHbTgov1d_x1)^2; ...                                % 1
    numel(idxHbTgov1d_x1)*numel(idxHbTgov1d_x2); ...            % 2
    numel(idxHbTgov1d_x2)^2; ...                                % 3
    numel(idxHbTgov1d_Pmech)*numel(idxHbTgov1d_x1); ...         % 4
    numel(idxHbTgov1d_Pmech)^2];                                % 5
assert(isequal(nBlkVals.',repmat(4,1,5)))
% Running start offset of each block inside a JvalsLin column
offBlk = [0; cumsum(nBlkVals)];
nValsLin = offBlk(end);
assert(nValsLin==5*4)
assert(nValsLin==20)
JvalsLin = zeros(nValsLin,n);

% State 1 (backpropgation for x2)
% x1 + T3*dx1/dt - x2 - T2*dx2/dt = 0
% (Block 1) d/dx1
blk = stampReImPages(1 + 1j*wk*T3);
JvalsLin(offBlk(1)+(1:nBlkVals(1)),:) = reshape(blk,nBlkVals(1),n);
% (Block 2) d/dx2
blk = stampReImPages(-1 - 1j*wk*T2);
JvalsLin(offBlk(2)+(1:nBlkVals(2)),:) = reshape(blk,nBlkVals(2),n);

% State 2 (backpropagation for pref)
% pref - (wr/ws-1) - R*x2 - R*T1*dx2/dt = 0
% (Block 3) d/dx2
blk = stampReImPages(-R - 1j*wk*R.*T1);   % R and T1 are pages, so ".*"
JvalsLin(offBlk(3)+(1:nBlkVals(3)),:) = reshape(blk,nBlkVals(3),n);

% Pmech
% x1 - (wr/ws-1)*Dt - pmech = 0
% (Block 4) d/dx1. Constant, so it is the same column for every governor.
blk = stampReImPages(1);
JvalsLin(offBlk(4)+(1:nBlkVals(4)),:) = repmat(blk(:),1,n);
% (Block 5) d/dPmech. Constant.
blk = stampReImPages(-1);
JvalsLin(offBlk(5)+(1:nBlkVals(5)),:) = repmat(blk(:),1,n);

% ------------------- Coupling to the GENROU speed Wr ---------------------
% The -(wr/ws) terms sit in the GENROU Wr columns, so they leave the
% governor block and are returned separately. ws is the NOMINAL speed
% (the legacy's constant obj.w), not the solved fundamental w0, which is why
% these two blocks are harmonic-independent.
% Rows 3:4 of the governor, from "-(wr/ws-1)" in state 2
dFx2_dWr = repmat(reshape(stampReImPages(-1/w0nom),[],1),1,n);
% Rows 5:6 of the governor, from "-(wr/ws-1)*Dt" in the Pmech equation
dFpmech_dWr = reshape(stampReImPages(-Dt/w0nom),4,n);

% ------------------------- Linear residual -------------------------------
% RESIDUAL SIGN. The legacy assembles Y*x = I. The caller here works with
% J*x - F, i.e. J*x = F, and J is the legacy's Ygov, so F is the legacy's
% Igov with NO sign change. The terms with no state attached are the "+pref"
% and the "+1" inside -(wr/ws-1); both move to the right-hand side with a
% flipped sign, giving -Pref-1 on state 2 and -Dt on the Pmech equation -
% exactly the legacy's "Igov(3) = -Pref - 1*N" and "Igov(5) = -Dt*N".
% They land on the real rows (3 and 5) and only at k=0 because both are pure
% DC quantities.
% N SCALING. The HB state vector holds DFT coefficients scaled by N (the
% iDFT divides by N), so every DC constant must carry that N. Pref arrives
% already multiplied (the caller builds gov_Pref0 = pref*N), but the literal
% "1" of -(wr/ws-1) and the bare Dt are raw time-domain constants, so the N
% is applied to them here. That is why this kernel takes N and the exciter
% kernel does not.
F = zeros(nEqns,n);
if isDc
    F(3,:) = -reshape(Pref,1,n) - 1*N;  % Only real idx
    F(5,:) = -reshape(Dt,1,n) * N;      % Only real idx
end

end


function Y_ReIm = stampReImPages(Y)
%STAMPREIMPAGES  Page-aware stampReIm, i.e. [real(Y) -imag(Y); imag(Y) real(Y)]
%   written with the real and imaginary rows/cols interleaved. Identical to
%   stampReIm.m, extended with a trailing page dimension so one call covers
%   every governor.

[nr,nc,n] = size(Y);
Y_ReIm = zeros(2*nr,2*nc,n);

Y_ReIm(1:2:end,1:2:end,:) =  real(Y);
Y_ReIm(1:2:end,2:2:end,:) = -imag(Y);
Y_ReIm(2:2:end,1:2:end,:) =  imag(Y);
Y_ReIm(2:2:end,2:2:end,:) =  real(Y);

end
