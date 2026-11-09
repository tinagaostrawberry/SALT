function [JvalsTm,FTm] = calcHbTgov1dMatrices_FD( ...
    ... % Time domain previous guesses (pages: [nVar x N x numTgov1d])
    pmechTd,TmTd,wrTd, ...
    ... % System constants
    w0nom,genP, ...
    ... % Fourier transform things
    K,GAMMA,invGAMMA)
%CALCHBTGOV1DMATRICES_FD  Vectorised TGOV1D nonlinear harmonic-balance stamp.
%
%   The single nonlinear TGOV1D equation (mechanical torque) couples every
%   harmonic to every other one, so each Jacobian block is the dense DFT
%   sandwich (GAMMA .* dF/dx(t)) * invGAMMA over all K+1 retained
%   harmonics. Devices ride the third dimension so the whole fleet is built
%   in one pass.
%
%   JvalsTm is [12*(K+1)^2 x numTgov1d] for the three blocks written into
%   the GENROU Tm rows 33:34 (w.r.t. gov Pmech, gen Wr, gen Tm), each
%   flattened column-major. FTm is [2*(K+1) x numTgov1d].
%
%   Ported from calcGovMatrices_nonlinear_FD.m (legacy harmonic balance),
%   whose separate stampOfGovMatrices_nonlinear_FD.m placement step is now
%   the caller's job: this kernel only emits the dense blocks and residual.

% VARIABLE LAYOUT (fixed): one TGOV1D block per harmonic is 6 rows and one
% GENROU block is 34, every quantity a (real,imag) pair. Within a device the
% row index is k*nEqns + localIndex, i.e. harmonic-outer / variable-inner -
% the ordering get_idxTermRow_ReIm.m builds in the legacy, and the same
% ordering hbDft produces when its output is read column-major.
% INDICIES NOTATATION:
% Var name          Indices (Re,Im)
% x1                1:2
% x2                3:4
% Pmech             5:6
idxHbTgov1d_Pmech = 5:6;
% GENROU rows/cols this kernel writes into (placed by the caller)
idxHbGenrou_Wr = 23:24;
idxHbGenrou_Tm = 33:34;

Kt = K + 1;

% Governors run along dim 3; time runs along dim 2.
n = size(pmechTd,3);
assert(all(cellfun(@(x)size(x,3)==n,{TmTd,wrTd})))

% ------------------------- Tm in generator -------------------------------
% Check electrical and mechanical angles are equal
% theta = (p/2)*thetaM (Section 3.1.2 in Kundar's textbook), which means
% wr = (p/2)*wM, since wr = dtheta/dt. If p = 2, then wr = (2/2)*wM = wM.
assert(genP==2)
% When wr = wM, we have:
% pmech - (wr/w)*Tm = 0
% (Deriv)
% d/dPmech has no time dependence, so it is the identity rather than a
% sandwich; repmat gives it the device dimension the other two blocks
% already have so all three reshape the same way.
dFgovTm_dPmech = repmat(eye(2*Kt),1,1,n);
dFgovTm_dWr    = dftSandwich(GAMMA,invGAMMA,-TmTd/w0nom);
dFgovTm_dTm    = dftSandwich(GAMMA,invGAMMA,-wrTd/w0nom);
% (Val)
FTm = hbDftPages(GAMMA,K, ...
    pmechTd - (wrTd/w0nom).*TmTd);

% Blocks into the GENROU Tm rows, in triplet-map order:
%   1. gen 33:34 x gov 5:6    (dFgovTm_dPmech)
%   2. gen 33:34 x gen 23:24  (dFgovTm_dWr)
%   3. gen 33:34 x gen 33:34  (dFgovTm_dTm)
% All three are square because every row and column set is a single
% (Re,Im) pair spread over the K+1 harmonics.
nValsTmBlock = numel(idxHbGenrou_Tm)*Kt*numel(idxHbTgov1d_Pmech)*Kt;
assert(nValsTmBlock == (2*Kt)^2)
assert(isequal(numel(idxHbGenrou_Wr),numel(idxHbGenrou_Tm)))
JvalsTm = [ ...
    reshape(dFgovTm_dPmech,nValsTmBlock,n); ...
    reshape(dFgovTm_dWr,nValsTmBlock,n); ...
    reshape(dFgovTm_dTm,nValsTmBlock,n)];
assert(size(JvalsTm,1) == 3*nValsTmBlock)

end

function M = dftSandwich(GAMMA,invGAMMA,dTd)
% Dense cross-harmonic Jacobian block of a nonlinearity whose time-domain
% derivative row is dTd ([1 x N x nDev]).
% The broadcast .* IS the diagonal multiplication GAMMA*diag(dTd); keeping
% it as a broadcast avoids ever forming the N-by-N diagonal.
% pagemtimes broadcasts the 2-D invGAMMA across every device page.
M = pagemtimes(GAMMA .* dTd, invGAMMA);
end

function F = hbDftPages(GAMMA,K,xTd)
% Residual of a scalar-per-device time-domain signal xTd ([1 x N x nDev]),
% returned as [2*(K+1) x nDev].
nDev = size(xTd,3);
% hbDft treats dim 1 as the list of signals to transform and the transform
% is independent per signal, so put the devices on dim 1 and do the whole
% fleet in one call.
X = hbDft(GAMMA,K,permute(xTd,[3 2 1]));   % [2*nDev x (K+1)]
% hbDft interleaves (Re,Im) down the rows per signal and puts harmonics
% along the columns; each device's column must instead run
% harmonic-outer / (Re,Im)-inner, hence split (Re,Im) from device, swap
% harmonic ahead of device, then collapse.
F = reshape(permute(reshape(X,2,nDev,K+1),[1 3 2]),2*(K+1),nDev);
end
