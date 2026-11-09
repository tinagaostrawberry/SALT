function [JvalsLin,F,dParkV_dVphA,dFdeltaTheta_dW0] = stampOfHbGenrouMatrices_FD( ...
    ... % Current guess (pages: [nState x 1 x numGenrou])
    i,e,efd,wr,deltaTheta,IphABC,Tm, ...
    ... % Harmonic and frequency
    k,w0, ...
    ... % Parameters (pages: [... x numGenrou])
    Rmat,Lmat,H,KD, ...
    ... % Constants
    w0nom,N, busToGenV,busToGenI)
%STAMPOFHBGENROUMATRICES_FD  Linear harmonic-balance GENROU stamp, one harmonic.
%
%   Harmonic-diagonal part of the GENROU equations, evaluated once per
%   harmonic k. Every input carries the machines along the third dimension,
%   so the whole fleet is stamped in one pass (the legacy looped one machine
%   at a time). JvalsLin is [nValsLin x numGenrou], one column per machine,
%   with the dense blocks flattened column-major and concatenated in the
%   fixed order documented below. F is the linear residual [34 x numGenrou]
%   for this harmonic. dParkV_dVphA and dFdeltaTheta_dW0 are the couplings
%   out of the device block: into the terminal bus voltage at the SAME
%   harmonic, and into the system-frequency column.
%
%   Ported from stampOfGenMatrices_linear_FD.m (legacy harmonic balance).
%   NOTE: the legacy also overwrote the imaginary DC rows of the device with
%   an identity to keep the matrix full rank. That is an assembly-level
%   overwrite, not a stamp, so it stays with the caller.

% VARIABLE LAYOUT (fixed): one GENROU block per harmonic is 34 rows, every
% quantity a (real,imag) pair. Within a device the row index is
% k*34 + localIndex, i.e. harmonic-outer / variable-inner.
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
idxHbGenrou_Efd        = 21:22;
idxHbGenrou_Wr         = 23:24;
idxHbGenrou_DeltaTheta = 25:26;
idxHbGenrou_IphABC     = 27:32;
idxHbGenrou_Tm         = 33:34;
nEqns = 34;
% Stator + field rows take the e = [ed eq e0 | efd] columns together
idxIdqf = idxHbGenrou_Idq0(1:8);
idxEdqf = [idxHbGenrou_Edq0, idxHbGenrou_Efd];

% Machines run along dim 3. The linear part has no state dependence, so the
% current guess is only used to pin down the fleet size and to check that
% every page count agrees.
n = size(wr,3);
assert(all(cellfun(@(x)size(x,3)==n, ...
    {i,e,efd,deltaTheta,IphABC,Tm,Rmat,Lmat,H,KD})))

% Harmonic frequency. k is an integer harmonic index, so the DC test is
% exact here and replaces the legacy's abs(wk)<1e-12 tolerance.
wk = k*w0;
isDc = (k==0);

% Block sizes, in the order they are concatenated into JvalsLin:
%  # | rows  | cols  | value
%  1 | 1:8   | 15:22 | stampReIm(-eye(4))                        ->  64
%  2 | 1:14  | 1:14  | stampReIm(-R + 1j*wk*L_PU)                -> 196
%  3 | 23:24 | 23:24 | stampReIm(1j*wk*2*H/w0nom + KD/w0nom)     ->   4
%  4 | 23:24 | 33:34 | stampReIm(-1)                             ->   4
%  5 | 25:26 | 23:24 | stampReIm(1)                              ->   4
%  6 | 25:26 | 25:26 | stampReIm(-1j*wk)                         ->   4
%  7 | 27:32 | 27:32 | stampReIm(-eye(3))*busToGenI           ->  36
% NOTE: the port spec wrote this total as 308, but 64+196+4+4+4+4+36 = 312.
% The seven blocks above are authoritative, so the count asserted here is
% 312; the 308 was an arithmetic slip in the spec.
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = [ ...
    numel(idxIdqf)*numel(idxEdqf); ...                      % 1
    numel(idxHbGenrou_Idq0)^2; ...                          % 2
    numel(idxHbGenrou_Wr)^2; ...                            % 3
    numel(idxHbGenrou_Wr)*numel(idxHbGenrou_Tm); ...        % 4
    numel(idxHbGenrou_DeltaTheta)*numel(idxHbGenrou_Wr); ...% 5
    numel(idxHbGenrou_DeltaTheta)^2; ...                    % 6
    numel(idxHbGenrou_IphABC)^2];                           % 7
assert(isequal(nBlkVals.',[64 196 4 4 4 4 36]))
% Running start offset of each block inside a JvalsLin column
offBlk = [0; cumsum(nBlkVals)];
nValsLin = offBlk(end);
assert(nValsLin==312)
JvalsLin = zeros(nValsLin,n);

% Branch eqns from Kundar's power syst book (Eqns 3.120-3.134):
% e = -[R]*i + [L]*(d/dt)i + [Lsp]*i*wr
%                            |_________|
%                           NONLIN, IGNORE
% s.t.: i = [id iq i0 | ifd i1d | i1q i2q]^T
%       e = [ed eq e0 | efd 0   | 0     0]^T
%       [Lsp]*i*wr = [-phiq*wr phid*wr | 0 ]
%       [L]*i = phi
%
% The first three rows are the dq0 stator windings, whose inductance is
% referred to nominal speed; the rotor windings are already in those units.
L_PU = Lmat./[w0nom;w0nom;w0nom;1;1;1;1];
% (Block 1) e. Constant, so it is the same column for every machine.
blk_dEdqf = stampReImPages(-eye(4));
JvalsLin(offBlk(1)+(1:nBlkVals(1)),:) = repmat(blk_dEdqf(:),1,n);
% (Block 2) i
blk_dI = stampReImPages(-Rmat + 1j*wk*L_PU);
JvalsLin(offBlk(2)+(1:nBlkVals(2)),:) = reshape(blk_dI,nBlkVals(2),n);

% Park transformation, Kundar's power txtbook (3.60):
% vabc = [T(w0*t + deltaTheta)^-1]*vdq0
%        |____________________________|
%                NONLIN, IGNORE
% Couples rows 15:20 to the terminal bus voltage at the SAME harmonic, so it
% leaves the device block and is returned separately. busToGenV is a
% per-machine scalar, hence ".*" rather than the legacy's "*"; the repmat
% forces the page dimension so a fleet-wide scalar still reshapes per page.
dParkV_dVphA_pages = repmat(stampReImPages(-eye(3)),1,1,n).*busToGenV;
dParkV_dVphA = reshape(dParkV_dVphA_pages,36,n);

% Equation for efd
% Depends on exc control, so leave row 21:22 of Jacobian empty for now!

% Swing eqns from Kundar's power syst book (Eqn 3.204):
% (d/dt)(deltaWr) = (1/(2*H))*(Tm - Te - KD*deltaWr)
% s.t.: Te = -phiq*id + phid*iq = i^T*([-phiq phid | 0 ])
%          = i^T*Lsp*i
% =>
% (2*H/w)*(d/dt)(wr) + (KD/w)*wr - Tm + i^T*Lsp*i - (KD/w)*w = 0
%                                       |_______|   |_______|
%                                         NONLIN,     CONST
%                                         IGNORE
% (Block 3) d(swing)/d(wr)
blk_dWr = stampReImPages(1j*wk*2*H/w0nom + KD/w0nom);
JvalsLin(offBlk(3)+(1:nBlkVals(3)),:) = reshape(blk_dWr,nBlkVals(3),n);
% (Block 4) d(swing)/d(Tm). Constant.
blk_dTm = stampReImPages(-1);
JvalsLin(offBlk(4)+(1:nBlkVals(4)),:) = repmat(blk_dTm(:),1,n);

% deltaWr = (1/w)*(wr - w0) = (1/w)*(d(theta)/dt - w0)
% => d(theta)/dt = wr
%
%
% Given that theta := w0*t + deltaTheta
% => deltaTheta = theta - w0*t
% => d(deltaTheta)/dt = d(theta)/dt - w0
% => d(deltaTheta)/dt = wr - w0
%
% NOTE: Should NOT replace this equation with deltaTheta=0 because
% deltaTheta at DC is just theta_DC!! Will be used accordingly in
% nonlinear equation
% (Block 5) d(deltaTheta eqn)/d(wr). Constant, stamped at every harmonic.
blk_dDeltaTheta_dWr = stampReImPages(1);
JvalsLin(offBlk(5)+(1:nBlkVals(5)),:) = repmat(blk_dDeltaTheta_dWr(:),1,n);
% (Block 6) d(deltaTheta eqn)/d(deltaTheta). The legacy skipped this stamp
% at DC; stampReIm(-1j*wk) is exactly zero there, so one branchless
% expression reproduces both legacy branches.
blk_dDeltaTheta = stampReImPages(-1j*wk);
JvalsLin(offBlk(6)+(1:nBlkVals(6)),:) = repmat(blk_dDeltaTheta(:),1,n);
% Coupling into the system-frequency (w0) column
if isDc
    % wr defined by DFT, xN forces w0 to be too
    dFdeltaTheta_dW0 = repmat([-1;0]*N,1,n);
else
    % No w0 at higher harmonics by defn of fundamental
    dFdeltaTheta_dW0 = zeros(2,n);
end

% ...and do same Park for current
% iabc = [T(w0*t + deltaTheta)^-1]*idq0
%        |____________________________|
%              NONLIN, IGNORE
% (Block 7) busToGenI is a per-machine scalar, hence ".*"; the repmat
% forces the page dimension so a fleet-wide scalar still reshapes per page.
blk_dIphABC = repmat(stampReImPages(-eye(3)),1,1,n).*busToGenI;
JvalsLin(offBlk(7)+(1:nBlkVals(7)),:) = reshape(blk_dIphABC,nBlkVals(7),n);

% Equation for TorqueM
% Depends on gov control, so leave row 33 to 34 of Jacobian empty for now!

% ------------------------- Linear residual -------------------------------
% RESIDUAL SIGN. The legacy assembles Y*x = I ("Simulation (Y*x_v = I, so I
% is negative)", EmtHb.m). The caller here works with J*x - F, i.e. J*x = F,
% and J is the legacy's Ygen, so F is the legacy's Igen with NO sign change.
% The only linear term with no state attached is the swing equation's
% -(KD/w)*w = -KD, which moves to the right-hand side as +KD. It is a pure
% DC quantity, and the HB state vector holds DFT coefficients scaled by N
% (the iDFT divides by N), so its DC coefficient is +KD*N. It lands on the
% real row of wr (row 23) and only at harmonic k=0 - exactly the legacy's
% "Igen(23) = KD * N" under "if isDc".
F = zeros(nEqns,n);
if isDc
    F(23,:) = N*reshape(KD,1,n);
end

end


function Y_ReIm = stampReImPages(Y)
%STAMPREIMPAGES  Page-aware stampReIm, i.e. [real(Y) -imag(Y); imag(Y) real(Y)]
%   written with the real and imaginary rows/cols interleaved. Identical to
%   stampReIm.m, extended with a trailing page dimension so one call covers
%   every machine.

[nr,nc,n] = size(Y);
Y_ReIm = zeros(2*nr,2*nc,n);

Y_ReIm(1:2:end,1:2:end,:) =  real(Y);
Y_ReIm(1:2:end,2:2:end,:) = -imag(Y);
Y_ReIm(2:2:end,1:2:end,:) =  imag(Y);
Y_ReIm(2:2:end,2:2:end,:) =  real(Y);

end
