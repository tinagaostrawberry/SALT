function [JvalsLin,F,dParkV_dVabc_lin,dFdw0,dFdw0Bare] = ...
    stampOfHbLoadCompMatrices_FD( ...
    ... % Current guess (pages: [nState x 1 x numLoadComp])
    vs,iphABC,wslip,deltaThetaSys,is,ir,wr,PQ,Vm,izip, ...
    ... % Harmonic and frequency
    k,w0, ...
    ... % Parameters (pages: [... x numLoadComp])
    R,L_PU,L_Te,T0,m,Tau, PQz,PQi,PQp, ...
    ... % Constants
    w0nom,N,H, busToLoadCompV,loadCompToBusI)
%STAMPOFHBLOADCOMPMATRICES_FD  Linear harmonic-balance composite-load stamp.
%
%   Harmonic-diagonal part of the composite load (induction motor + ZIP),
%   evaluated once per harmonic k. Every input carries the loads along the
%   third dimension, so the whole fleet is stamped in one pass. JvalsLin is
%   [nValsLin x numLoadComp], one column per load, with the dense blocks
%   flattened column-major and concatenated in the fixed order documented
%   below. F is the linear residual [34 x numLoadComp] for this harmonic.
%   dParkV_dVabc_lin is deliberately empty (see below); dFdw0 and dFdw0Bare
%   are the two halves of the coupling into the system-frequency column -
%   see "Coupling into the system-frequency column" for why the caller must
%   keep them apart when it builds the right-hand side.
%
%   Generalised from stampOfCompositeLoadMatrices_SS.m by the rule in the
%   HB class header: take SALT's algebraic form verbatim, then add back each
%   dropped derivative as (j*k*w0) acting on the harmonic coefficient,
%   keeping whatever constant already stood in front of d/dt. At k = 0 every
%   restored term vanishes and this reduces to SALT exactly.
%
%   TIME BASE. Physical seconds throughout, as the TD model
%   (stampOfCompositeLoadMatrices_TD_flat.m) discretises it: the stator /
%   rotor inductance enters as L_PU = L/w0nom, the swing inertia as
%   2*H/w0nom, and the ZIP magnitude filter as Tau in seconds. So every
%   restored derivative is j*k*w0 times those same constants.
%
%   TWO STRUCTURAL CHANGES FROM SALT, both load-bearing:
%
%     (a) The Park transforms are REAL, not short-circuited. SALT sets
%         theta_sys = 0 and collapses vdq_s = T(theta_sys)*vabc to
%         vdq_s = vphA, which only holds at one harmonic. HB uses the
%         genuine transform with theta_sys(t) = w0*t + deltaThetaSys(t),
%         which makes both Park rows a DFT-sandwich nonlinearity exactly
%         like the GENROU's. They therefore live in
%         calcHbLoadCompMatrices_FD, not here.
%
%     (b) SALT's FIphBC = iphBC - phShiftBC*iphA is GONE. That is a
%         positive-sequence condition (two complex constraints); the model
%         actually assumes only i_a + i_b + i_c = 0. The terminal current is
%         instead written abc FROM dq through the 3x2 inverse Park, whose
%         columns [cos(theta_abc), -sin(theta_abc)] each already sum to zero,
%         so the zero-sum holds by construction and needs no equation of its
%         own - see the zero-sum / abc-from-dq paragraph in HBINITINDICES.
%         The voltage goes the other way, dq FROM abc through the 2x3
%         forward Park.
%
%   See also CALCHBLOADCOMPMATRICES_FD, STAMPOFCOMPOSITELOADMATRICES_SS,
%   STAMPOFHBGENROUMATRICES_FD.

% VARIABLE LAYOUT (fixed): one composite-load block per harmonic is 34
% rows, every quantity a (real,imag) pair, so a dq pair costs 4 slots and a
% scalar 2. Within a device the row index is k*34 + localIndex, i.e.
% harmonic-outer / variable-inner, and a dq pair is stored
% (d_re, d_im, q_re, q_im).
% Var name          Indices (Re,Im)
% Vs (vd_s,vq_s)     1:4
% IphABC (a,b,c)     5:10      BUS per-unit
% WSlip             11:12
% DeltaThetaSys     13:14
% Is (id_s,iq_s)    15:18
% Ir (id_r,iq_r)    19:22
% Wr                23:24
% PQ (P,Q)          25:28
% Vm                29:30
% Izip (id,iq)      31:34
idxHbLoadComp_Vs            =  1:4;
idxHbLoadComp_IphABC        =  5:10;
idxHbLoadComp_WSlip         = 11:12;
idxHbLoadComp_DeltaThetaSys = 13:14;
idxHbLoadComp_Is            = 15:18;
idxHbLoadComp_Ir            = 19:22;
idxHbLoadComp_Wr            = 23:24;
idxHbLoadComp_PQ            = 25:28;
idxHbLoadComp_Vm            = 29:30;
idxHbLoadComp_Izip          = 31:34;
nEqns = 34;
% The stator/rotor equation is one 4x4 complex block over isr = [is; ir], so
% its rows and columns are the two index vectors taken together.
idxIsr = [idxHbLoadComp_Is, idxHbLoadComp_Ir];
% ...and only the STATOR half of that block sees the terminal voltage,
% because vsr = [vs; 0; 0].
idxSrStator = idxHbLoadComp_Is;

% Loads run along dim 3. The linear part has no state dependence, so the
% current guess is only used to pin down the fleet size and to check that
% every page count agrees.
n = size(wr,3);
assert(all(cellfun(@(x)size(x,3)==n, ...
    {vs,iphABC,wslip,deltaThetaSys,is,ir,PQ,Vm,izip, ...
     R,L_PU,L_Te,T0,m,Tau,PQz,PQi,PQp,H, ...
     busToLoadCompV,loadCompToBusI})))
% Every state page is a column of (re,im) slots in the layout above
assert(size(vs,1)==numel(idxHbLoadComp_Vs))
assert(size(iphABC,1)==numel(idxHbLoadComp_IphABC))
assert(size(wslip,1)==numel(idxHbLoadComp_WSlip))
assert(size(deltaThetaSys,1)==numel(idxHbLoadComp_DeltaThetaSys))
assert(size(is,1)==numel(idxHbLoadComp_Is))
assert(size(ir,1)==numel(idxHbLoadComp_Ir))
assert(size(wr,1)==numel(idxHbLoadComp_Wr))
assert(size(PQ,1)==numel(idxHbLoadComp_PQ))
assert(size(Vm,1)==numel(idxHbLoadComp_Vm))
assert(size(izip,1)==numel(idxHbLoadComp_Izip))
% N is the DFT length, one scalar for the whole system, not a per-load page
assert(isscalar(N))
% R, L_PU and L_Te are the 4x4 (d_s,q_s,d_r,q_r) matrices SALT builds
assert(size(R,1)==4 && size(R,2)==4)
assert(size(L_PU,1)==4 && size(L_PU,2)==4)
assert(size(L_Te,1)==4 && size(L_Te,2)==4)

% Harmonic frequency. k is an integer harmonic index, so the DC test is
% exact.
wk = k*w0;
isDc = (k==0);

% Block sizes, in the order they are concatenated into JvalsLin:
%  # | rows  | cols  | value                                      | count
%  1 | 1:4   | 1:4   | -eye(4)                                    |  16
%  2 | 5:10  | 5:10  | -eye(6)                                    |  36
%  3 | 11:12 | 11:12 | -eye(2)                                    |   4
%  4 | 11:12 | 23:24 | -eye(2)                                    |   4
%  5 | 13:14 | 13:14 | -eye(2)                                    |   4
%  6 | 15:18 | 1:4   | -eye(4)  (only the stator rows see vs)     |  16
%  7 | 15:22 | 15:22 | stampReIm(R + 1j*wk*L_PU + wMatSys*L_PU)   |  64
%  8 | 23:24 | 23:24 | stampReIm(-(2*H/w0nom)*1j*wk)              |   4
%  9 | 25:28 | 25:28 | -eye(4)                                    |  16
% 10 | 29:30 | 29:30 | stampReIm(-(1 + 1j*wk*Tau))                |   4
% 11 | 31:34 | 31:34 | -eye(4)                                    |  16
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = [ ...
    numel(idxHbLoadComp_Vs)^2; ...                                  %  1
    numel(idxHbLoadComp_IphABC)^2; ...                              %  2
    numel(idxHbLoadComp_WSlip)^2; ...                               %  3
    numel(idxHbLoadComp_WSlip)*numel(idxHbLoadComp_Wr); ...         %  4
    numel(idxHbLoadComp_DeltaThetaSys)^2; ...                       %  5
    numel(idxSrStator)*numel(idxHbLoadComp_Vs); ...                 %  6
    numel(idxIsr)^2; ...                                            %  7
    numel(idxHbLoadComp_Wr)^2; ...                                  %  8
    numel(idxHbLoadComp_PQ)^2; ...                                  %  9
    numel(idxHbLoadComp_Vm)^2; ...                                  % 10
    numel(idxHbLoadComp_Izip)^2];                                   % 11
assert(isequal(nBlkVals.',[16 36 4 4 4 16 64 4 16 4 16]))
% Running start offset of each block inside a JvalsLin column
offBlk = [0; cumsum(nBlkVals)];
nValsLin = offBlk(end);
assert(nValsLin==184)
JvalsLin = zeros(nValsLin,n);

% ------------------------- Induction motor -------------------------------

% Park V (rows 1:4)
% SALT:  vdq_s = T(theta_sys)*vabc  =>  vdq_s = vphA
% HB:    -vs + T2(thetaSys)*v_abc*busToLoadCompV = 0
% Only the -vs half is linear; T2(thetaSys)*v_abc is a DFT sandwich and is
% stamped by calcHbLoadCompMatrices_FD (its Jacobian w.r.t. deltaThetaSys
% and w.r.t. the bus voltages, and its residual).
% (Block 1) d(ParkV)/d(vs). Constant, so the same column for every load.
blk_dParkV_dVs = -eye(numel(idxHbLoadComp_Vs));
JvalsLin(offBlk(1)+(1:nBlkVals(1)),:) = repmat(blk_dParkV_dVs(:),1,n);

% Park I (rows 5:10)
% SALT:  (idq_s+idq_ZIP) = T(theta_sys)*iabc  =>  (idq_s+idq_ZIP) = iphA
% HB:    -i_abc + invT2(thetaSys)*(is + izip)*loadCompToBusI = 0
% Written abc FROM dq, so the +/-120 degree balance equations disappear:
% each column of invT2 = [cos(theta_abc), -sin(theta_abc)] sums to zero over
% the three phases, hence i_a+i_b+i_c = 0 by construction. Only the -i_abc
% half is linear. NOTE the direction flip of the base change relative to
% SALT: iphABC is in BUS per-unit and the dq currents are in LOAD per-unit,
% so the inverse Park carries loadCompToBusI where SALT's forward Park
% carried busToLoadCompI.
% (Block 2) d(ParkI)/d(iphABC). Constant.
blk_dParkI_dIphABC = -eye(numel(idxHbLoadComp_IphABC));
JvalsLin(offBlk(2)+(1:nBlkVals(2)),:) = repmat(blk_dParkI_dIphABC(:),1,n);

% Slip frequency (rows 11:12)
% SALT:  wSlip = wSys - wr
% HB:    -wslip + w0*delta_k0 - wr = 0
% wSys is the system fundamental, a SINGLE SCALAR unknown, so it is not
% harmonic-expanded: it couples only at k = 0, exactly as the GENROU's
% d(deltaTheta)/dt = wr - w0 does ("no w0 at higher harmonics by definition
% of the fundamental"). wr IS harmonic-expanded, so it is stamped at every
% harmonic. There is no d/dt in this equation, so nothing is restored.
% (Block 3) d(Fslip)/d(wslip) and (Block 4) d(Fslip)/d(wr). Both constant.
blk_dFslip_dWslip = -eye(numel(idxHbLoadComp_WSlip));
JvalsLin(offBlk(3)+(1:nBlkVals(3)),:) = repmat(blk_dFslip_dWslip(:),1,n);
blk_dFslip_dWr = -eye(numel(idxHbLoadComp_WSlip));
JvalsLin(offBlk(4)+(1:nBlkVals(4)),:) = repmat(blk_dFslip_dWr(:),1,n);

% Frame angle (rows 13:14)
% SALT:  d(thetasys)/dt = wSys => thetasys = wSys*t => thetasys = 0 @ t=0
% HB:    -deltaThetaSys = 0
% HB writes theta_sys = w0*t + deltaThetaSys, and the system frame rotates
% at exactly the fundamental (wSys == w0), so d(deltaThetaSys)/dt = 0 and
% the DC value is SALT's "theta_sys = 0 at t = 0". Both statements are the
% single equation deltaThetaSys_k = 0 at every k. Unlike the GENROU's
% deltaTheta, this one really is pinned: the GENROU angle is a machine state
% driven by wr - w0, this one is the definition of the frame.
% (Block 5) d(Fsys)/d(deltaThetaSys). Constant.
blk_dFsys = -eye(numel(idxHbLoadComp_DeltaThetaSys));
JvalsLin(offBlk(5)+(1:nBlkVals(5)),:) = repmat(blk_dFsys(:),1,n);

% Stator/rotor equations (rows 15:22)
% SALT:
% |vs| = [R]*|idq_s| + (d/dt)*[L]/w0*|idq_s| + [wMat]/w0*[L]*|idq_s|
% |0 |       |idq_r|                 |idq_r|                 |idq_r|
% =>
% |vs| = [R]*|idq_s| + [wMat]/w0*[L]*|idq_s|
% |0 |       |idq_r|                 |idq_r|
%
% s.t.
% R = |Rs_I2x2      |          L = |Lss_I2x2 Lm_I2x2 |
%     |      Rr_I2x2|              |Lm_I2x2  Lrr_I2x2|
%
% HB restores the dropped derivative as j*k*w0 acting on the coefficient,
% with the SAME constant L_PU = L/w0nom in front of it (physical seconds):
%   -vsr + R*isr + j*k*w0*L_PU*isr + wMat*L_PU*isr = 0
%
% wMat SPLITS IN TWO, and the split is the whole reason this model needs two
% kernels:
%   stator block (rows 1:2) uses wSys, which is DC-ONLY and equals the
%       scalar unknown w0. So wMatSys*L_PU*isr is LINEAR in isr with
%       coefficient w0 and is harmonic-diagonal - it belongs here, plus a
%       d/dw0 column.
%   rotor block (rows 3:4) uses wslip, which IS harmonic-expanded, so
%       wMatSlip*L_PU*isr is a product of two unknowns and goes through the
%       DFT sandwich in calcHbLoadCompMatrices_FD.
wMatSys = zeros(4,4,n);
wMatSys(1,2,:) = -w0;  wMatSys(2,1,:) =  w0;
% (Block 6) d/d(vs) of -vsr. vsr = [vs; 0; 0], so only the stator rows see
% the terminal voltage and the block is 4x4 real, not 8x4.
blk_dSr_dVs = -eye(numel(idxSrStator));
JvalsLin(offBlk(6)+(1:nBlkVals(6)),:) = repmat(blk_dSr_dVs(:),1,n);
% (Block 7) d/d(isr). stampReIm of the complex 4x4 turns it into the 8x8
% real block with (re,im) interleaved, which is exactly rows/cols 15:22.
blk_dSr_dIsr = stampReImPages(R + 1j*wk*L_PU + pagemtimes(wMatSys,L_PU));
JvalsLin(offBlk(7)+(1:nBlkVals(7)),:) = reshape(blk_dSr_dIsr,nBlkVals(7),n);

% Swing equation (rows 23:24)
% SALT:
% (2*H/w0)*d(wr)/dt = ((isr^T)*L_Te*isr) - T0*(wr/wsys)^m
% =>
% 0 = ((isr^T)*L_Te*isr) - T0*(wr/wsys)^m
% HB restores the derivative with the same 2*H/w0nom in front (physical
% seconds, as the TD model discretises it):
%   -(2*H/w0nom)*(j*k*w0)*wr + isr'*L_Te*isr - T0*(wr/wSys)^m = 0
% The sign is SALT's with the derivative moved to the same side as the rest
% of the residual, hence the leading minus. Both remaining terms are
% nonlinear (a quadratic form and a power of wr) and are stamped by
% calcHbLoadCompMatrices_FD; wSys there is again the scalar w0, so the load
% torque also needs a d/dw0 column (SALT's dFswing_dwSys).
% (Block 8) d/d(wr) of the restored inertia term.
blk_dSwing_dWr = stampReImPages(-(2*H/w0nom)*1j*wk);
JvalsLin(offBlk(8)+(1:nBlkVals(8)),:) = reshape(blk_dSwing_dWr,nBlkVals(8),n);

% ---------------------------- ZIP model ----------------------------------

% ZIP power (rows 25:28)
% SALT:  PQ = PQz*Vm^2 + PQi*Vm + PQp
% HB:    -PQ + PQz*Vm^2 + PQi*Vm + PQp = 0        (no d/dt to restore)
% Only the -PQ half is linear. The Vm^2 and Vm terms are nonlinear in Vm, so
% they go through the sandwich; PQp is a pure constant and rides along in
% the nonlinear kernel's residual value, where the DFT puts it in the DC bin
% with the N scaling the state vector uses. That is why PQp is unused here
% (it is still taken as an argument so both kernels share one parameter
% list, and so its page count is checked above).
% (Block 9) d(FZip)/d(PQ). Constant.
blk_dZip_dPQ = -eye(numel(idxHbLoadComp_PQ));
JvalsLin(offBlk(9)+(1:nBlkVals(9)),:) = repmat(blk_dZip_dPQ(:),1,n);

% ZIP voltage magnitude (rows 29:30)
% SALT:  Vm = (vd_s^2+vq_s^2)^(1/2)
% The SS kernel has NO filter because d/dt = 0 makes Vm = |vs|; the TD model
% carries d(Vm)/dt = (1/Tau)*(|vs| - Vm). Restoring that derivative and
% multiplying through by Tau gives
%   -(1 + j*k*w0*Tau)*Vm + |vs| = 0
% so the ZIP sees an attenuated voltage at every harmonic and behaves as an
% impedance there. At k = 0 this collapses back to SALT's -Vm + |vs|.
% |vs| is nonlinear and is stamped by calcHbLoadCompMatrices_FD.
% (Block 10) d(FVm)/d(Vm).
blk_dVm = stampReImPages(-(1 + 1j*wk*Tau));
JvalsLin(offBlk(10)+(1:nBlkVals(10)),:) = reshape(blk_dVm,nBlkVals(10),n);

% ZIP current (rows 31:34)
% SALT:  idZip = (P*vd + Q*vq)/Vm^2
%        iqZip = (P*vq - Q*vd)/Vm^2
% HB:    -izip + [vs, J*vs]*PQ / Vm^2 = 0,  J = [0 1; -1 0]
% (no d/dt to restore). Only the -izip half is linear.
% (Block 11) d(FIZip)/d(izip). Constant.
blk_dIzip = -eye(numel(idxHbLoadComp_Izip));
JvalsLin(offBlk(11)+(1:nBlkVals(11)),:) = repmat(blk_dIzip(:),1,n);

% ----------------- Park V coupling into the bus voltage -------------------
% RETURNED EMPTY ON PURPOSE. For the GENROU the terminal-voltage coupling is
% an identity and therefore harmonic-diagonal and linear, so its kernel can
% hand it back here. The composite load's Park V multiplies the bus voltage
% by T2(thetaSys), which is a DFT sandwich dense across harmonics, so the
% whole coupling belongs in calcHbLoadCompMatrices_FD (its JvalsVabc). The
% output slot is kept so the signature mirrors the GENROU's and the caller
% can be written against one shape.
dParkV_dVabc_lin = [];

% ------------------------- Linear residual -------------------------------
% RESIDUAL SIGN AND WHY F IS IDENTICALLY ZERO.
% The caller works with J*x - F, i.e. J*x = F. For a purely linear equation
% that collapses to MINUS the equation's state-free constant: the GENROU
% settles it, where the swing's -(KD/w)*w = -KD constant comes back as
% F = +KD*N (HB.stampGenrou, "for a LINEAR equation the Newton right-hand
% side J*x - F(x) collapses to the equation's CONSTANT term"). The N is
% there because the HB state vector holds UNNORMALIZED DFT coefficients, so
% a DC physical constant c appears as c*N.
%
% SIGN OF THE SLIP EQUATION'S w0 TERM, derived (not guessed) from SALT:
% stampOfCompositeLoadMatrices_SS has Fslip = -wSlip + wSys - wr with
% dFslip_dwSys = +1, and the TD model agrees (Fslip = -wsl + wS - wr). So
% wSys - here the scalar unknown w0 - enters the slip residual with a PLUS
% sign, at k = 0 only. That is the only term in these 34 equations with no
% device state attached.
%
% BUT it is not a constant: w0 is the solver's frequency unknown and it has
% its own Jacobian column, returned below as dFdw0 (+1 on the slip row at
% k = 0, the same +1 SALT stamps). HB.stampGenrou states explicitly that
% "the linear J*x must NOT be added here - and neither must the terminal-
% voltage column or the w0 column, which are linear too and cancel in the
% same way". The +w0 therefore cancels between J*x and the residual exactly
% like every other linear term, and contributes nothing to F.
%
% Net result: there is no state-free constant anywhere in the 34 linear
% rows (PQp is the only other candidate and it is carried in the nonlinear
% kernel's residual, where the DFT applies the N scaling), so F is all
% zeros. The slot is kept because the caller adds it unconditionally. If a
% future caller chooses NOT to stamp the w0 column for the slip row, the
% value it would need there instead is F(11,:) = -w0*N at k = 0 - the minus
% because F is MINUS the constant, as the GENROU's -KD -> +KD shows.
F = zeros(nEqns,n);

% ----------------- Coupling into the system-frequency column -------------
% d/dw0 of every row, at this harmonic. Four groups contribute.
%
% WHERE THE N GOES. w0 is a plain scalar unknown, while every other unknown
% in x_v is a DFT coefficient, and at DC a coefficient equals N times the
% time-domain value. So:
%   a BARE ADDITIVE w0 sits alongside DFT coefficients and needs *N - this
%       is exactly why the GENROU's angle equation carries [-1;0]*N;
%   a d/dw0 of a term like 1j*k*w0*L_PU*isr does NOT get N, because the N is
%       already carried inside isr.
% Audited row by row below: only the slip row (i) has a bare additive w0.
%
% THE COLUMN IS RETURNED IN TWO PIECES, and the split is not cosmetic - it
% is what keeps the Newton right-hand side exact. The caller's convention
% (HB.stampGenrou) is that J*x is added to the right-hand side only for
% NONLINEAR terms, because for a linear term J*x and F cancel. A bare
% additive w0 is linear in w0, so its column cancels and must NOT be added;
% every other entry here is the w0 half of a term that is BILINEAR in w0 and
% a state (1j*k*w0*L_PU*isr) or outright nonlinear in w0 (the load torque),
% for which J*x does NOT cancel and the column's contribution dFdw0*w0 has
% to go into the right-hand side. Hence:
%   dFdw0Bare  group (i) only      - stamp in J, keep OUT of the RHS
%   dFdw0      groups (ii)-(iv)    - stamp in J, and add dFdw0*w0 to the RHS
% The Jacobian's w0 column is their sum.
%
% The GENROU needs no such split because its only w0 column IS group (i):
% it deliberately omits the bilinear w0 columns at k > 0, which leaves J*x
% equal to F there and the right-hand side at zero.
dFdw0     = zeros(nEqns,n);
dFdw0Bare = zeros(nEqns,n);

% (i) Slip row: BARE ADDITIVE w0, so +N, at k = 0 only. The +1 is straight
% from SALT's dFslip_dwSys = +1 (Fslip = -wSlip + wSys - wr, and the TD model
% agrees); the N is the coefficient convention above, the same "wr defined by
% DFT, xN forces w0 to be too" as the GENROU's dFdeltaTheta_dW0 = [-1;0]*N.
% w0 is a real DC scalar, so it lands on the REAL row of wslip (row 11) and
% never on the imaginary one.
if isDc
    dFdw0Bare(idxHbLoadComp_WSlip(1),:) = N;
end

% (ii) Stator/rotor rows: d/dw0 of (1j*k*w0*L_PU*isr + wMatSys*L_PU*isr).
% The restored derivative contributes 1j*k*L_PU*isr on ALL four equations;
% the speed-voltage term contributes only on the stator pair and is exactly
% SALT's dFs_dwSys = [-phisr(2); phisr(1)].
% NO N: both terms multiply isr, whose coefficients already carry it.
isrC   = toComplexPages([is; ir]);          % [4 x 1 x n] complex coefficients
phisrC = pagemtimes(L_PU,isrC);             % flux linkage, same units as SALT
dSr_dw0 = 1j*k*phisrC;
dSr_dw0(1,:,:) = dSr_dw0(1,:,:) - phisrC(2,:,:);
dSr_dw0(2,:,:) = dSr_dw0(2,:,:) + phisrC(1,:,:);
dFdw0(idxIsr,:) = reshape(toReImPages(dSr_dw0),numel(idxIsr),n);

% (iii) Swing row: d/dw0 of the restored inertia term
% -(2*H/w0nom)*1j*k*w0*wr. Algebraic in the coefficient, so it belongs here.
% NO N: it multiplies wr, whose coefficients already carry it.
%
% THE LOAD TORQUE'S w0 COLUMN IS NOT HERE, DELIBERATELY. SALT has
% dFswing_dwSys = -T0*wr^m*(-m)*wSys^(-m-1), and writing that same algebra
% on the harmonic COEFFICIENT of wr - which is all this kernel is handed -
% is wrong by N^(m-1), because a DC coefficient is N times the time value
% and the m-th power multiplies that error m times over while the single
% 1/w0 divides it out once. Pointwise in time the derivative is just
% -(m/w0) times the load torque itself, so the exact column is -(m/w0)
% times its DFT, which only calcHbLoadCompMatrices_FD can form. It returns
% it as its own dFdw0 and the caller adds the two columns. The error was not
% cosmetic: at m = 2, N = 3 it is a factor of 3 on this entry, which cost
% Newton its quadratic convergence and, since this column also feeds the
% right-hand side, moved the converged answer.
wrC = toComplexPages(wr);                   % [1 x 1 x n]
dSwing_dw0 = -(2*H/w0nom).*(1j*k).*wrC;
dFdw0(idxHbLoadComp_Wr,:) = ...
    reshape(toReImPages(dSwing_dw0),numel(idxHbLoadComp_Wr),n);

% (iv) ZIP voltage row: d/dw0 of -(1j*k*w0*Tau)*Vm. The "1" and the |vs| in
% that equation carry no w0.
% NO N: the term multiplies Vm, whose coefficients already carry it.
VmC = toComplexPages(Vm);                   % [1 x 1 x n]
dVm_dw0 = -1j*k*Tau.*VmC;
dFdw0(idxHbLoadComp_Vm,:) = ...
    reshape(toReImPages(dVm_dw0),numel(idxHbLoadComp_Vm),n);

% Rows 19:22 get nothing from (ii)'s speed-voltage half and rows 1:10,
% 13:14, 25:28 and 31:34 nothing at all: the rotor speed term is driven by
% wslip (a state, handled in the nonlinear kernel), and the Park, ZIP power
% and ZIP current equations contain no w0.
%
% BARE-CONSTANT AUDIT over all 34 rows, for the same reason N appears above.
% The only state-free constant anywhere in this equation set besides the
% slip row's w0 is the ZIP power row's PQp. It is carried in the nonlinear
% kernel's residual value, where the DFT of a constant lands in the DC bin
% already multiplied by N - so it is correctly scaled without being touched
% here. Every other term is a product of unknowns, so nothing else needs N.

% UNUSED ARGUMENTS, all deliberate and all page-count-checked above:
%   L_Te             the electric torque is a quadratic form, so it is
%                    nonlinear and belongs to the other kernel.
%   T0, m            the load torque and its only Jacobian column both need
%                    the TIME-domain wr, so they live in the other kernel
%                    too (see the swing row's note in the w0 section). Both
%                    are still taken here so one parameter list serves both
%                    kernels and their page counts are checked above - and m
%                    is read by nothing else, which is the point.
%   PQp              a pure constant, carried in the other kernel's
%                    residual so the DFT applies the N scaling (see above).
%   busToLoadCompV   the Park V coupling multiplies by T2(thetaSys), so the
%                    base change travels with it into the other kernel.
%   loadCompToBusI   likewise, the base change rides on the inverse Park;
%                    the linear half of the Park I row is a bare -eye(6),
%                    unlike the GENROU's, which carries busToGenI because
%                    its abc current needs no transform here.
%   iphABC, deltaThetaSys, PQ, izip, and the state values generally, which
%                    the linear part does not depend on at all.

end


function Y_ReIm = stampReImPages(Y)
%STAMPREIMPAGES  Page-aware stampReIm, i.e. [real(Y) -imag(Y); imag(Y) real(Y)]
%   written with the real and imaginary rows/cols interleaved. Identical to
%   stampReIm.m, extended with a trailing page dimension so one call covers
%   every load.

[nr,nc,n] = size(Y);
Y_ReIm = zeros(2*nr,2*nc,n);

Y_ReIm(1:2:end,1:2:end,:) =  real(Y);
Y_ReIm(1:2:end,2:2:end,:) = -imag(Y);
Y_ReIm(2:2:end,1:2:end,:) =  imag(Y);
Y_ReIm(2:2:end,2:2:end,:) =  real(Y);

end


function Xc = toComplexPages(x)
%TOCOMPLEXPAGES  Fold an interleaved (re,im) state page into complex values.
%   x is [2*nVar x 1 x nDev] with the real part of each quantity on the odd
%   rows and the imaginary part on the even ones - the layout the whole
%   solver uses - so the stride-2 slices ARE the two parts.

Xc = x(1:2:end,:,:) + 1j*x(2:2:end,:,:);

end


function x = toReImPages(Xc)
%TOREIMPAGES  Inverse of TOCOMPLEXPAGES: unfold complex values back to the
%   interleaved (re,im) row layout, so the result drops straight onto the
%   device's rows.

[nr,nc,n] = size(Xc);
x = zeros(2*nr,nc,n);
x(1:2:end,:,:) = real(Xc);
x(2:2:end,:,:) = imag(Xc);

end
