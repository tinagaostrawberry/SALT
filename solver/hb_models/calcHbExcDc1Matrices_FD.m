function [JvalsNlin,Fnlin,JvalsEfd,FEfd] = calcHbExcDc1Matrices_FD( ...
    ... % Time domain previous guesses (pages: [nVar x N x numExcDc1])
    x1Td,VtTd,vabcTd,efdTd,wrTd, ...
    ... % System constants
    w0nom,excDc1_genToExc,busToGenV, ...
    ... % Fourier transform things
    K,GAMMA,invGAMMA)
%CALCHBEXCDC1MATRICES_FD  Vectorised EXDC1 nonlinear harmonic-balance stamp.
%
%   The two nonlinear EXDC1 equations couple every harmonic to every other
%   one, so each Jacobian block is the dense DFT sandwich
%   (GAMMA .* dF/dx(t)) * invGAMMA over all K+1 retained harmonics. Devices
%   ride the third dimension so the whole fleet is built in one pass.
%
%   JvalsNlin is [16*(K+1)^2 x numExcDc1]: the terminal-voltage row's own
%   block (exc rows 11:12 x exc cols 11:12) followed by its bus-voltage
%   block (exc rows 11:12 x the six Vabc (Re,Im) bus columns), each
%   flattened column-major. Fnlin is [12*(K+1) x numExcDc1], the exciter's
%   own nonlinear residual (FexcVt on rows 11:12, zero elsewhere).
%
%   JvalsEfd is [12*(K+1)^2 x numExcDc1] for the three blocks the exciter
%   writes into the GENROU Efd rows 21:22 (w.r.t. exc x1, gen Efd, gen Wr),
%   and FEfd is [2*(K+1) x numExcDc1], that equation's residual.
%
%   Ported from calcExcMatrices_nonlinear_FD.m (legacy harmonic balance),
%   whose separate stampOfExcMatrices_nonlinear_FD.m placement step is now
%   the caller's job: this kernel only emits the dense blocks and residuals.
%
%   UNCERTAINTY FLAG: excDc1_genToExc and busToGenV carry the roles the legacy
%   called exc_genToExc and busToGenV respectively, i.e. excDc1_genToExc must be
%   gen_Lad/gen_Rfd (it multiplies the GENROU field voltage to put it in
%   exciter PU) and busToGenV must be the bus-PU -> gen-PU voltage
%   factor, NOT its reciprocal. The legacy names read the other way round
%   and the voltage factor is 1 in both the legacy and the SS kernels, so
%   its direction is not observable from the sources.

% VARIABLE LAYOUT (fixed): one EXDC1 block per harmonic is 12 rows and one
% GENROU block is 34, every quantity a (real,imag) pair. Within a device the
% row index is k*nEqns + localIndex, i.e. harmonic-outer / variable-inner -
% the ordering get_idxTermRow_ReIm.m builds in the legacy, and the same
% ordering hbDft produces when its output is read column-major.
% INDICIES NOTATATION:
% Var name          Indices (Re,Im)
% x1                1:2
% x2                3:4
% x3                5:6
% x4                7:8
% x5                9:10
% Vt                11:12
idxHbExcDc1_x1 =  1:2;
idxHbExcDc1_Vt = 11:12;
nExcEqns = 12;
% GENROU rows this kernel writes into (placed by the caller)
idxHbGenrou_Efd = 21:22;
idxHbGenrou_Wr  = 23:24;
% Bus voltage block: 3 phases, each a (Re,Im) pair, per harmonic
nVabcVars = 3;

Kt = K + 1;

% Exciters run along dim 3; time runs along dim 2.
n = size(x1Td,3);
assert(all(cellfun(@(x)size(x,3)==n, ...
    {VtTd,vabcTd,efdTd,wrTd,excDc1_genToExc,busToGenV})))
assert(size(vabcTd,1)==nVabcVars)

% ------------------------- Efd in generator ------------------------------
% Efd normalized by speed
% efd - x1*wr_pu = 0
% (Deriv)
dFEfd_dX1     = dftSandwich(GAMMA,invGAMMA,-wrTd/w0nom);
% d/defd has no time dependence, so it is a scaled identity rather than a
% sandwich; repmat gives it the device dimension the other two blocks
% already have so all three reshape the same way.
dFexcEfd_defd = repmat(eye(2*Kt),1,1,n).*excDc1_genToExc;
dFexcEfd_dwr  = dftSandwich(GAMMA,invGAMMA,-x1Td/w0nom);
% (Val)
FEfd = hbDftPages(GAMMA,K, ...
    efdTd.*excDc1_genToExc - x1Td.*wrTd/w0nom);

% Blocks into the GENROU Efd rows, in triplet-map order:
%   1. gen 21:22 x exc 1:2    (dFEfd_dX1)
%   2. gen 21:22 x gen 21:22  (dFexcEfd_defd)
%   3. gen 21:22 x gen 23:24  (dFexcEfd_dwr)
% All three are square because every row and column set is a single
% (Re,Im) pair spread over the K+1 harmonics.
nValsEfdBlock = numel(idxHbGenrou_Efd)*Kt*numel(idxHbExcDc1_x1)*Kt;
assert(nValsEfdBlock == (2*Kt)^2)
assert(isequal(numel(idxHbGenrou_Wr),numel(idxHbGenrou_Efd)))
JvalsEfd = [ ...
    reshape(dFEfd_dX1,nValsEfdBlock,n); ...
    reshape(dFexcEfd_defd,nValsEfdBlock,n); ...
    reshape(dFexcEfd_dwr,nValsEfdBlock,n)];
assert(size(JvalsEfd,1) == 3*nValsEfdBlock)

% ------------------------ Terminal voltage vt ----------------------------
% -Vt + sqrt(2/3)*(va^2 + vb^2 + vc^2)^(1/2) = 0
% s.t. va,vb,vc are gen PU
% (Const)
vtabcTd_genPU = busToGenV.*vabcTd;
% Sum over the three phases only: dim 1 must be named because the signal
% now carries time on dim 2 and devices on dim 3.
vtabc_sumOfSq_genPU = sum(vtabcTd_genPU.*vtabcTd_genPU,1);
% Column positions of phase jdxVar inside the dense Vabc block. Columns run
% harmonic-outer / variable-inner, so the harmonic stride is the 6 columns
% of one bus voltage triple and the phase contributes a (Re,Im) pair offset.
% This is the legacy's get_idxTermRow_ReIm(jdxVar,3,K) written out.
idxVtCol_ReIm = @(jdxVar)reshape( ...
    (0:K)*2*nVabcVars + ((2*jdxVar-1):(2*jdxVar)).',[],1);
% (Deriv)
dFexcVt_dVt = repmat(-eye(2*Kt),1,1,n);
dFexcVt_dVabc = zeros(2*Kt,2*nVabcVars*Kt,n);
for jdxVar = 1:nVabcVars
    dFexcVt_dVabc(:,idxVtCol_ReIm(jdxVar),:) = ...
        dftSandwich(GAMMA,invGAMMA, ...
        sqrt(2/3)*( vtabc_sumOfSq_genPU.^(-1/2) ).*vtabcTd_genPU(jdxVar,:,:)) ...
        .*busToGenV; % This conversion is for when multiplying deriv by vabc
end
% (Val)
FexcVt = hbDftPages(GAMMA,K, ...
    -VtTd + sqrt(2/3)*(vtabc_sumOfSq_genPU.^(1/2)));

% Blocks into the exciter's own Vt rows, in triplet-map order:
%   1. exc 11:12 x exc 11:12  (dFexcVt_dVt)
%   2. exc 11:12 x bus Vabc   (dFexcVt_dVabc)
% The bus columns use the caller's own spacing, so only the dense block is
% emitted here and the caller places it.
nBlkVals = [ ...
    numel(idxHbExcDc1_Vt)*Kt*numel(idxHbExcDc1_Vt)*Kt; ...  % 1
    numel(idxHbExcDc1_Vt)*Kt*2*nVabcVars*Kt];               % 2
nValsNlin = sum(nBlkVals);
assert(nValsNlin == 16*Kt^2)
JvalsNlin = [ ...
    reshape(dFexcVt_dVt,nBlkVals(1),n); ...
    reshape(dFexcVt_dVabc,nBlkVals(2),n)];
assert(size(JvalsNlin,1) == nValsNlin)

% The exciter's own residual is zero everywhere except the Vt equation.
% Rows run harmonic-outer / variable-inner with stride nExcEqns, which is
% why the harmonic offset multiplies nExcEqns and not the pair width.
Fnlin = zeros(nExcEqns*Kt,n);
Fnlin(reshape((0:K)*nExcEqns + idxHbExcDc1_Vt.',[],1),:) = FexcVt;

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
