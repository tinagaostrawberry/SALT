function [JvalsNlin,Fnlin] = calcHbGenrouMatrices_FD( ...
    ... % Time-domain previous guess (pages: [nVar x N x numGenrou])
    iTd,vTd,wrTd,deltaThetaTd, ...
    ... % Frequency-domain previous guess, used for J*x
    w0, ...
    ... % Parameters (pages: [... x numGenrou])
    Lspmat, w0nom, ...
    ... % Transform
    K,N,GAMMA,invGAMMA,tSamp)
%CALCHBGENROUMATRICES_FD  Nonlinear harmonic-balance GENROU kernel.
%
%   Dense-across-harmonics part of the GENROU equations: each derivative is
%   formed in the time domain and squeezed back into harmonic coefficients
%   by the DFT sandwich (GAMMA .* dRow) * invGAMMA. Every input carries the
%   machines along the third dimension, so the whole fleet is evaluated in
%   one pass (the legacy looped one machine at a time). JvalsNlin is
%   [nValsNlin x numGenrou] with the dense blocks flattened column-major and
%   concatenated in the fixed order documented below; Fnlin is the nonlinear
%   residual [34*(K+1) x numGenrou].
%
%   Ported from calcGenMatrices_nonlinear_FD.m and the RHS half of
%   stampOfGenMatrices_nonlinear_FD.m (legacy harmonic balance).

% VARIABLE LAYOUT (fixed): one GENROU block per harmonic is 34 rows, every
% quantity a (real,imag) pair, so a global row is k*34 + localIndex, i.e.
% harmonic-outer / variable-inner. The legacy's get_idxTermRow_ReIm builds
% exactly this ordering inside each dense block (stride numFreqEqns*numVars
% per harmonic, (re,im) innermost), so the legacy blocks need no reordering.
% Var #         Var name        Indices (Re,Im)
% 1:7           Idq0            1:14
% 8:10          Edq0            15:20
% 11            Efd             21:22
% 12            Wr              23:24
% 13            deltaTheta      25:26
% 14:16         IphABC          27:32
% 17            Tm              33:34
idxHbGenrou_Idq0       =  1:14;
idxHbGenrou_Edq0       = 15:20;
idxHbGenrou_Wr         = 23:24;
idxHbGenrou_DeltaTheta = 25:26;
idxHbGenrou_IphABC     = 27:32;
nEqns = 34;
% The current Park transform only sees the dq0 stator currents
idxHbGenrou_Idq0_abcOnly = idxHbGenrou_Idq0(1:6);

% K is # harmonics, NOT including zeros, so total is +1
K_tot = K + 1;
% Machines run along dim 3
n = size(wrTd,3);
% Time runs along dim 2. getInvParkMatrices_concatTheta asserted that theta
% was a row vector; the page-based equivalent is that the scalar states
% carry one row and tSamp is a column of samples.
assert(size(wrTd,1)==1 && size(deltaThetaTd,1)==1 && size(tSamp,2)==1)
assert(size(iTd,2)==N && size(vTd,2)==N && size(tSamp,1)==N)
% Each DFT sandwich must come out [2*K_tot x 2*K_tot] to land on a harmonic
% block, which holds only when GAMMA keeps exactly DC..K (N = 2*K+1).
assert(size(GAMMA,1)==2*K_tot && size(GAMMA,2)==N)
assert(size(invGAMMA,1)==N && size(invGAMMA,2)==2*K_tot)

% Block layout, in the order they are concatenated into JvalsNlin. Each is a
% dense (numel(rows)*K_tot) x (numel(cols)*K_tot) array.
%  # | rows  | cols  | legacy name
%  1 | 1:14  | 1:14  | dSpeedVolt_dI           -> 196*K_tot^2
%  2 | 1:14  | 23:24 | dSpeedVolt_dWr          ->  28*K_tot^2
%  3 | 15:20 | 15:20 | dInvParkV_dV            ->  36*K_tot^2
%  4 | 15:20 | 25:26 | dInvParkV_dDeltaTheta   ->  12*K_tot^2
%  5 | 23:24 | 1:14  | dTe_dI                  ->  28*K_tot^2
%  6 | 27:32 | 1:6   | dInvParkI_dI            ->  36*K_tot^2
%  7 | 27:32 | 25:26 | dInvParkI_dDeltaTheta   ->  12*K_tot^2
nBlkVals = K_tot^2 * [ ...
    numel(idxHbGenrou_Idq0)*numel(idxHbGenrou_Idq0); ...
    numel(idxHbGenrou_Idq0)*numel(idxHbGenrou_Wr); ...
    numel(idxHbGenrou_Edq0)*numel(idxHbGenrou_Edq0); ...
    numel(idxHbGenrou_Edq0)*numel(idxHbGenrou_DeltaTheta); ...
    numel(idxHbGenrou_Wr)*numel(idxHbGenrou_Idq0); ...
    numel(idxHbGenrou_IphABC)*numel(idxHbGenrou_Idq0_abcOnly); ...
    numel(idxHbGenrou_IphABC)*numel(idxHbGenrou_DeltaTheta)];
assert(isequal(nBlkVals.',K_tot^2*[196 28 36 12 28 36 12]))
% Running start offset of each block inside a JvalsNlin column
offBlk = [0; cumsum(nBlkVals)];
nValsNlin = offBlk(end);
assert(nValsNlin==348*K_tot^2)
JvalsNlin = zeros(nValsNlin,n);

% Number of (re,im) pairs addressed by each block, i.e. the "numTerms" the
% legacy passed to get_idxTermRow_ReIm
nVarsI  = numel(idxHbGenrou_Idq0)/2;             % 7
nVarsE  = numel(idxHbGenrou_Edq0)/2;             % 3
nVarsAbc = numel(idxHbGenrou_IphABC)/2;          % 3

% -------------------- Voltage speed --------------------------
% speedVolt = [Lsp]*i*wr
% (Const)
Lsp_PU = Lspmat/w0nom;
% (Deriv)
dSpeedVolt_dI  = zeros(2*nVarsI*K_tot,2*nVarsI*K_tot,n);
dSpeedVolt_dWr = zeros(2*nVarsI*K_tot,2*K_tot,n);
for idxVar = 1:2 % Speed terms only apply to ed and eq
    for jdxVar = 1:nVarsI
        dSpeedVolt_dI(idxTermRow_ReIm(idxVar,nVarsI,K), ...
            idxTermRow_ReIm(jdxVar,nVarsI,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            Lsp_PU(idxVar,jdxVar,:).*wrTd);
    end
    % Lsp_PU row times the whole current vector, per time sample and machine
    dSpeedVolt_dWr(idxTermRow_ReIm(idxVar,nVarsI,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA, ...
        pagemtimes(Lsp_PU(idxVar,:,:),iTd));
end
% (Val)
SpeedVolt = pagemtimes(Lsp_PU,iTd).*wrTd;
% ------------------ Inverse Park (V,I) -----------------------
% NOTE: Need to double-check!! Since w0 is fundamental frequency, there
% should be a derivative term w.r.t. w0. However, its role in the Park
% transform is purely enforcing periodicty (e.g. cos(w0*t + deltaTheta)) so
% in the frequency domain, the change of w0 doesn't affect the fourier
% coefficients of the outcome of the Park transform...
%
% invParkV = [invPark(w0*t + theta_DC + deltaTheta)]*v
% invParkI = [invPark(w0*t + theta_DC + deltaTheta)]*i
% (Const)
% tSamp is transposed so time lines up with dim 2 of the page arrays
thetaTd = w0*(tSamp.') + deltaThetaTd;
% Inlined from getInvParkMatrices_concatTheta: that helper concatenated the
% three dq0 columns along the (row-vector) time axis, which cannot carry a
% page dimension, so the entries are kept as explicit [3 x N x n] rows here.
thetaTd_phAbc = thetaTd + [0;-2*pi/3;2*pi/3];
cosTd = cos(thetaTd_phAbc);
sinTd = sin(thetaTd_phAbc);
% invT  row a = [ cos(theta_a) -sin(theta_a) 1 ]
% JinvT row a = [-sin(theta_a) -cos(theta_a) 0 ]
invTcols  = {cosTd, -sinTd, ones(nVarsAbc,N,n)};
% (Deriv + Val)
dInvParkV_dV          = zeros(2*nVarsE*K_tot,2*nVarsE*K_tot,n);
dInvParkV_dDeltaTheta = zeros(2*nVarsAbc*K_tot,2*K_tot,n);
dInvParkI_dDeltaTheta = zeros(2*nVarsAbc*K_tot,2*K_tot,n);
for idxVar = 1:nVarsAbc
    for jdxVar = 1:nVarsE
        dInvParkV_dV(idxTermRow_ReIm(idxVar,nVarsE,K), ...
            idxTermRow_ReIm(jdxVar,nVarsE,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            invTcols{jdxVar}(idxVar,:,:));
    end
    % d(invT*x)/d(deltaTheta) contracted over dq0 at each time sample
    JinvPark_V = -sinTd(idxVar,:,:).*vTd(1,:,:) ...
        - cosTd(idxVar,:,:).*vTd(2,:,:);
    JinvPark_I = -sinTd(idxVar,:,:).*iTd(1,:,:) ...
        - cosTd(idxVar,:,:).*iTd(2,:,:);
    dInvParkV_dDeltaTheta(idxTermRow_ReIm(idxVar,nVarsAbc,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA,JinvPark_V);
    dInvParkI_dDeltaTheta(idxTermRow_ReIm(idxVar,nVarsAbc,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA,JinvPark_I);
end
% The inverse Park matrix is the same for voltage and current, so the legacy
% reuses the voltage derivative verbatim
dInvParkI_dI = dInvParkV_dV;
% (Val) invT row a dotted with the dq0 vector at each time sample
InvParkV = cosTd.*vTd(1,:,:) - sinTd.*vTd(2,:,:) + vTd(3,:,:);
InvParkI = cosTd.*iTd(1,:,:) - sinTd.*iTd(2,:,:) + iTd(3,:,:);
% ------------------- Electric Torque -------------------------
% Te = i^T*Lsp*i
% (Const)
% NOTE: kept as in the legacy - the torque uses the un-normalised Lsp while
% the speed voltage above uses Lsp/w0nom.
dte_di = pagemtimes(Lspmat + pagetranspose(Lspmat),iTd);
% (Deriv)
dTe_dI = zeros(2*K_tot,2*nVarsI*K_tot,n);
for jdxVar = 1:nVarsI
    dTe_dI(:,idxTermRow_ReIm(jdxVar,nVarsI,K),:) = ...
        dftSandwich(GAMMA,invGAMMA,dte_di(jdxVar,:,:));
end
% (Val) the quadratic form is per time sample, so contract dim 1 rather than
% calling pagemtimes with a time-indexed middle dimension
Te = sum(iTd.*pagemtimes(Lspmat,iTd),1);

% ------------------------ Jacobian values --------------------------------
JvalsNlin(offBlk(1)+(1:nBlkVals(1)),:) = ...
    reshape(dSpeedVolt_dI,nBlkVals(1),n);
JvalsNlin(offBlk(2)+(1:nBlkVals(2)),:) = ...
    reshape(dSpeedVolt_dWr,nBlkVals(2),n);
JvalsNlin(offBlk(3)+(1:nBlkVals(3)),:) = ...
    reshape(dInvParkV_dV,nBlkVals(3),n);
JvalsNlin(offBlk(4)+(1:nBlkVals(4)),:) = ...
    reshape(dInvParkV_dDeltaTheta,nBlkVals(4),n);
JvalsNlin(offBlk(5)+(1:nBlkVals(5)),:) = ...
    reshape(dTe_dI,nBlkVals(5),n);
JvalsNlin(offBlk(6)+(1:nBlkVals(6)),:) = ...
    reshape(dInvParkI_dI,nBlkVals(6),n);
JvalsNlin(offBlk(7)+(1:nBlkVals(7)),:) = ...
    reshape(dInvParkI_dDeltaTheta,nBlkVals(7),n);

% ----------------------- Nonlinear residual ------------------------------
% Fnlin holds the harmonic coefficients of the nonlinear terms themselves,
% g(x_prev), at the rows the terms belong to. The legacy folded these into
% Igen as (dg/dx)*x_prev - g(x_prev) (see stampOfGenMatrices_nonlinear_FD).
% With the caller's J*x - F convention, J*x = F, so the caller forms
% F = JvalsNlin*x_prev - Fnlin; the minus on g(x_prev) therefore lives in
% the caller and Fnlin is returned with the plain sign of the term.
Fnlin = zeros(nEqns*K_tot,n);
% Rows are gathered harmonic-outer / variable-inner, matching the layout
% that the DFT sandwich and hbDft both produce
Fnlin(hbRows(idxHbGenrou_Idq0,K,nEqns),:)   = dftPagesFlat(GAMMA,K,SpeedVolt);
Fnlin(hbRows(idxHbGenrou_Edq0,K,nEqns),:)   = dftPagesFlat(GAMMA,K,InvParkV);
Fnlin(hbRows(idxHbGenrou_Wr,K,nEqns),:)     = dftPagesFlat(GAMMA,K,Te);
Fnlin(hbRows(idxHbGenrou_IphABC,K,nEqns),:) = dftPagesFlat(GAMMA,K,InvParkI);

end


function S = dftSandwich(GAMMA,invGAMMA,dRow)
%DFTSANDWICH  (GAMMA .* dRow) * invGAMMA for every machine at once.
%   dRow is [1 x N x numGenrou]: one time-domain derivative row per machine.
%   The broadcast ".*" IS the diagonal multiplication by the derivative, so
%   it must never become an explicit diag. The result is
%   [2*(K+1) x 2*(K+1) x numGenrou], rows/cols ordered (re,im) within
%   ascending harmonic, which is the per-harmonic block ordering used here.

S = pagemtimes(GAMMA.*dRow,invGAMMA);

end


function rows = idxTermRow_ReIm(idxVar,numTerms,K)
%IDXTERMROW_REIM  Row indices of one variable inside a dense harmonic block.
%   Reproduces get_idxTermRow_ReIm with EmtHb.numFreqEqns fixed at 2: the
%   harmonic stride is 2*numTerms and (re,im) sit innermost, i.e.
%   harmonic-outer / variable-inner.

numFreqEqns = 2;
matSize = numFreqEqns*numTerms;
idxTermRow = (0:matSize:(matSize*K));
rows = reshape( ...
    idxTermRow + (idxVar-1)*numFreqEqns + (1:numFreqEqns).', ...
    [],1);

end


function rows = hbRows(locRows,K,nEqns)
%HBROWS  Device rows of a variable group across all harmonics.
%   Column-major reshape of (0:K)*nEqns + locRows(:) puts the harmonic
%   outermost and the variable innermost, matching the device row layout.

rows = reshape((0:K)*nEqns + locRows(:),[],1);

end


function Xflat = dftPagesFlat(GAMMA,K,xTd)
%DFTPAGESFLAT  Harmonic coefficients of a page-based time-domain signal.
%   xTd is [nVar x N x numGenrou]; the result is
%   [2*nVar*(K+1) x numGenrou], rows harmonic-outer / variable-inner with
%   (re,im) innermost.

[nVar,Nt,n] = size(xTd);
% GAMMA acts only along time, so the pages are folded into hbDft's variable
% dimension and the whole fleet goes through one 2-D call. Permuting to
% (var,page,time) keeps a machine's variables contiguous in the fold.
X2 = reshape(permute(xTd,[1 3 2]),nVar*n,Nt);
X_ReIm = hbDft(GAMMA,K,X2);             % [2*nVar*n x (K+1)]
% Undo the fold: within a column the rows run (re,im) inside variable,
% then machine, so splitting dim 1 recovers the page dimension
X3 = reshape(X_ReIm,2*nVar,n,K+1);
% Put harmonics ahead of machines, then flatten to harmonic-outer rows
Xflat = reshape(permute(X3,[1 3 2]),2*nVar*(K+1),n);

end
