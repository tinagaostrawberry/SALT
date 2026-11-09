function [JvalsLin,F,dFdw0,dFdw0Bare] = stampOfHbIbrGFMMatrices_FD( ...
    ... % Current guess (pages: [nState x 1 x numIbrGFM])
    IphABC,vg,vo,iff,io,vi,gamma,ifRef,x,wr,voRefDroopQ,pqavg, ...
    deltaTheta,voVRavg,voVR,voRef, ...
    ... % Harmonic and frequency
    k,w0, ...
    ... % Parameters (pages: [1 x 1 x numIbrGFM])
    Lf,rf,Cf,Rcap,Lc,rc,kc_i,kc_p,Gc,kv_i,kv_p,Gv,Mp,Mq, ...
    pset,qset,vo_q_set,wmeas_rad,w_VR_rad,Rv, ...
    ... % Transform
    ibrGFMToBusI, ...
    ... % Constants
    w0nom,N)
%STAMPOFHBIBRGFMMATRICES_FD  Linear HB stamp for a GFM inverter, one harmonic.
%
%   Harmonic-diagonal part of the grid-forming inverter equations, evaluated
%   once per harmonic k. Every input carries the inverters along the third
%   dimension, so the whole fleet is stamped in one pass. JvalsLin is
%   [nValsLin x numIbrGFM], one column per inverter, with the dense blocks
%   flattened column-major and concatenated in the fixed order documented
%   below. F is the linear residual [62 x numIbrGFM] for this harmonic
%   (CONSTANT TERMS ONLY - see "Linear residual"), and dFdw0 / dFdw0Bare are
%   the two halves of the coupling into the system-frequency column, each
%   [62 x numIbrGFM]; see "Coupling into the system-frequency column" for
%   why the caller must keep them apart when it builds the right-hand side.
%
%   Generalised from stampOfIbrGFMMatrices_SS.m by the rule in the HB class
%   header: take SALT's algebraic form verbatim, then add back each dropped
%   derivative as (j*k*w0) acting on the harmonic coefficient, keeping
%   whatever constant already stood in front of d/dt. At k = 0 every
%   restored term vanishes and this reduces to SALT exactly.
%
%   TIME BASE IS MIXED IN THIS DEVICE, and getting it wrong scales a
%   harmonic impedance by w0nom (~377) silently. The authority is the TD
%   model, stampOfIbrGFMMatrices_TD_flat.m, where the derivative is actually
%   discretised: dT is physical seconds, dTPU = dT*w0 is normalised time.
%   Verified block by block against that file:
%
%     LCL filter + coupling   dTPU  ("kC=-1/dTPU+...", "kP=-1/dTPU",
%     (io, vo, if rows)             "kF=-1/dTPU+...")
%                             -> restored derivative is j*k*w0/w0nom
%     inner current loop      dT    ("kG=kci/2+kcp/dT", "mdt=-1/dT")
%     inner voltage loop      dT    ("kV=kvi/2+kvp/dT", "mdt=-1/dT")
%     droop P/Q measurement   dT    ("Dp1=-1/dT+wm*(-1/2)")
%     frame angle             dT    ("Dt1=-1/dT")
%     virtual-resistor filter dT    ("Ra=-1/dT-wVR/2")
%                             -> restored derivative is j*k*w0
%
%   SIGN OF EVERY RESTORED DERIVATIVE also comes from the TD model, not from
%   guessing which side of SALT's residual it belongs on. The three passive
%   equations are written by SALT as "x' = RHS" and kept as "0 = RHS", so the
%   derivative comes back with a MINUS (TD: "FC1=-(iod-iodp)/dTPU+(1./Lc)*..."
%   and likewise FP1, FF1). The two integrator-error equations are the same
%   shape, so they are minus too (TD: "FKi1=-(gd-gdp)/dT+...",
%   "FVi1=-(xd-xdp)/dT+..."). The two inner-loop output equations carry their
%   derivative with a PLUS, because kc_p*d(gamma)/dt and kv_p*d(x)/dt sit on
%   the same side as the rest of the controller output (TD:
%   "FK1=...+kcp.*(gd-gdp)/dT+...", "FV1c=...+kvp.*(xd-xdp)/dT+...").
%
%   A PRE-EXISTING CONVENTION, REPRODUCED UNCHANGED: SALT divides both
%   measurement bandwidths by 2*pi (w_VR = w_VR_rad/(2*pi), and the same for
%   wmeas) in BOTH its steady-state and its time-domain kernel (TD:
%   "wm=wmR/(2*pi)", "wVR=wVRR/(2*pi)"), so it is a consistent convention
%   rather than a transcription slip - but a data value of 120 becomes a
%   ~19 rad/s filter. It matters more here than in SALT, because in HB these
%   two filters actually roll something off.
%
%   TWO NORMALISATIONS THAT LOOK INCONSISTENT AND ARE NOT. SALT carries the
%   virtual-resistor filter with its w_VR factor intact (FVRVoAvg =
%   w_VR*(Rv*io - vo_VR_avg)) but divides the power-measurement filter
%   through by wmeas (FDroopPavg = p - pavg). Restoring the derivative in
%   each case's OWN normalisation is therefore:
%       VR filter : -(j*k*w0)*vo_VR_avg      + w_VR*(Rv*io - vo_VR_avg) = 0
%       P filter  : -(j*k*w0/wmeas)*pavg     + p - pavg                 = 0
%   Both agree with the TD rows after dividing by the same factor SALT did.
%
%   THREE STRUCTURAL CHANGES FROM SALT:
%
%     (a) The Park transforms are REAL, not short-circuited. SALT uses a 2x2
%         invT that maps io_dq straight onto phase A, which only holds at one
%         harmonic. HB uses the genuine transforms with
%         theta(t) = w0*t + deltaTheta(t): the CURRENT goes dq -> abc through
%         the 3x2 inverse Park (3 complex = 6 rows), the VOLTAGE goes
%         abc -> dq through the 2x3 forward Park (2 complex = 4 rows). Both
%         are trig of an unknown angle times an unknown, so they are DFT
%         sandwiches and live in calcHbIbrGFMMatrices_FD, not here; only the
%         -i_abc and -vg_dq halves are linear.
%
%         PARK SIGN CONVENTION: getParkMatrices_2x3_PSCAD, which is what the
%         TD model inlines ("T11=2/3*c1 ... T21=2/3*s1"), uses
%         T = 2/3*[cos(theta_abc); sin(theta_abc)] with
%         theta_abc = theta + [0 -2*pi/3 2*pi/3], so its matching 3x2 inverse
%         is invT2 = [cos(theta_abc).', +sin(theta_abc).']. That is +sin,
%         which DIFFERS from the GENROU's getInvParkMatrices
%         ([cos, -sin, 1]). This device keeps the PSCAD convention, as SALT's
%         own 2x2 invT = [cos sin; sin -cos] already did.
%
%     (b) SALT's FIphBC = iphBC - phShiftBC*iphA is GONE. That is a
%         positive-sequence condition; the model actually assumes only
%         i_a + i_b + i_c = 0. Writing the terminal current abc FROM dq makes
%         the zero-sum hold by construction, because every column of
%         [cos(theta_abc); sin(theta_abc)] already sums to zero over the three
%         phases, so it needs no equation of its own.
%
%     (c) DeltaTheta, not Theta. SALT's FSpeed = -ws + wr came from
%         d(theta)/dt = wr. HB writes theta = w0*t + deltaTheta, so the
%         equation becomes d(deltaTheta)/dt = wr - w0, i.e.
%         -j*k*w0*deltaTheta_k + wr_k - w0*delta_k0 = 0, which at k = 0 pins
%         wr to w0 and recovers SALT. Signs and the w0 column mirror
%         stampOfHbGenrouMatrices_FD blocks 5 and 6 and its
%         dFdeltaTheta_dW0 exactly; the TD model agrees ("Dt1=-1/dT",
%         "Dt2=1/2", i.e. -theta' + wr).
%
%   TWO DERIVATIVES ARE RESTORED ON THE CAPACITOR ROW, NOT ONE. SALT's
%   capacitor equation is
%     vo_dq' = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq + Rcap*(if_dq' - io_dq')
%   and SALT drops BOTH vo_dq' AND the whole Rcap*(if' - io') term. Rcap is
%   the ESR in series with the filter capacitor - 0.005 pu in the 39-bus
%   data - so it is exactly the damping that decides the LCL resonance, and
%   dropping it would leave that resonance undamped. The TD model
%   discretises it in the SAME normalised time as the rest of the passive
%   block ("gI=(1./Cf)/2+Rcp/dTPU", "gO=(1./Cf)*(-1/2)+Rcp*(-1/dTPU)",
%   i.e. Rcap*(if'-io')/dT_PU), and with a PLUS
%   ("FP1=...+Rcp.*((ifd-ifdp)/dTPU-(iod-iodp)/dTPU)"), so it comes back as
%     + Rcap*(j*k*w0/w0nom)*(if_dq - io_dq)
%   It lands on the two capacitor blocks that already exist (rows 15:18 x
%   cols 15:18 and rows 15:18 x cols 19:22), so it is ADDED INTO blocks 7
%   and 8 rather than creating new ones: nValsLin is unchanged at 540 and
%   the block table below is unchanged. At k = 0 it vanishes with every
%   other restored derivative, so SALT is still recovered exactly.
%
%   THREE DEVIATIONS FROM THE PORT SPEC'S SIGNATURE:
%     - N is taken as a trailing constant (as the GENROU kernel does). The
%       linear residual has genuine state-free constants (the droop
%       setpoints), and the HB state vector holds UNNORMALIZED DFT
%       coefficients, so a DC physical constant c appears as c*N. Without N
%       those constants cannot be scaled.
%     - Rcap is taken next to Cf, for the term above.
%     - the spec's row table gives the three passive rows' restored
%       derivative as "+j*k*w0/w0nom*x"; the TD model says MINUS (see above),
%       and a one-branch inductor check settles it: 0 = -j*k*w0/w0nom*io +
%       (1/Lc)*(vo-vg) gives io = (vo-vg)/(j*k*w0*Lc_PU), the right
%       impedance, where the plus sign gives its negative.
%
%   See also CALCHBIBRGFMMATRICES_FD, STAMPOFIBRGFMMATRICES_SS,
%   STAMPOFHBGENROUMATRICES_FD, STAMPOFHBLOADCOMPMATRICES_FD.

% VARIABLE LAYOUT (fixed): one GFM-IBR block per harmonic is 62 rows, every
% quantity a (real,imag) pair, so a dq pair costs 4 slots and a scalar 2.
% Within a device the row index is k*62 + localIndex, i.e. harmonic-outer /
% variable-inner, and a dq pair is stored (d_re, d_im, q_re, q_im).
% Var name          Indices (Re,Im)
% IphABC (a,b,c)     1:6       BUS per-unit
% Vg (d,q)           7:10
% Vo                11:14
% If                15:18
% Io                19:22
% Vi                23:26
% Gamma             27:30
% IfRef             31:34
% X                 35:38
% Wr                39:40
% VoRefDroop_q      41:42
% PqAvg (pavg,qavg) 43:46
% DeltaTheta        47:48
% Vo_VR_avg         49:52
% Vo_VR             53:56
% VoRef             57:60
idxHbIbrGFM_IphABC       =  1:6;
idxHbIbrGFM_Vg           =  7:10;
idxHbIbrGFM_Vo           = 11:14;
idxHbIbrGFM_If           = 15:18;
idxHbIbrGFM_Io           = 19:22;
idxHbIbrGFM_Vi           = 23:26;
idxHbIbrGFM_Gamma        = 27:30;
idxHbIbrGFM_IfRef        = 31:34;
idxHbIbrGFM_X            = 35:38;
idxHbIbrGFM_Wr           = 39:40;
idxHbIbrGFM_VoRefDroop_q = 41:42;
idxHbIbrGFM_PqAvg        = 43:46;
idxHbIbrGFM_DeltaTheta   = 47:48;
idxHbIbrGFM_Vo_VR_avg    = 49:52;
idxHbIbrGFM_Vo_VR        = 53:56;
idxHbIbrGFM_VoRef        = 57:60;
idxHbIbrGFM_Io0          = 61:62;
nEqns = 62;
% The two halves of the PqAvg pair are addressed separately, because P and Q
% are two independent first-order filters.
idxHbIbrGFM_Pavg = idxHbIbrGFM_PqAvg(1:2);
idxHbIbrGFM_Qavg = idxHbIbrGFM_PqAvg(3:4);

% EQUATION LAYOUT (fixed, and NOT the same permutation as the variables -
% e.g. the io differential equation sits on rows 11:14 while the io variable
% occupies columns 19:22). Row order, with the SALT residual each one
% generalises:
% rows     equation                  SALT residual
%  1:6     Park I (abc from dq)      FParkI, with FIphBC removed
%  7:10    Park V (dq from abc)      FParkV
% 11:14    passive coupling          FPassiveCoup  - (j*k*w0/w0nom)*io
% 15:18    passive capacitor         FPassiveCap   - (j*k*w0/w0nom)*vo
% 19:22    passive filter            FPassiveFilt  - (j*k*w0/w0nom)*if
% 23:26    current controller        FibrCurrCont  + kc_p*(j*k*w0)*gamma
% 27:30    current integrator error  FibrCurrIntErr- (j*k*w0)*gamma
% 31:34    voltage controller        FVoltCont     + kv_p*(j*k*w0)*x
% 35:38    voltage integrator error  FVoltIntErr   - (j*k*w0)*x
% 39:40    P droop                   FDroopW
% 41:42    Q droop                   FDroopV
% 43:44    measured P                FDroopPavg    - (j*k*w0/wmeas)*pavg
% 45:46    measured Q                FDroopQavg    - (j*k*w0/wmeas)*qavg
% 47:48    frame angle               FSpeed, as d(deltaTheta)/dt = wr - w0
% 49:52    VR filter                 FVRVoAvg      - (j*k*w0)*vo_VR_avg
% 53:56    VR drop                   FVRVo
% 57:60    VR reference              FVRVoRef
% 61:62    zero sequence             v0 = (va+vb+vc)/3 = 0, the delta
%                                    connection; see CALCHBIBRGFMMATRICES_FD
rowIbrGFM_ParkI    =  1:6;
rowIbrGFM_ParkV    =  7:10;
rowIbrGFM_Coup     = 11:14;
rowIbrGFM_Cap      = 15:18;
rowIbrGFM_Filt     = 19:22;
rowIbrGFM_CurrCont = 23:26;
rowIbrGFM_CurrInt  = 27:30;
rowIbrGFM_VoltCont = 31:34;
rowIbrGFM_VoltInt  = 35:38;
rowIbrGFM_DroopW   = 39:40;
rowIbrGFM_DroopV   = 41:42;
rowIbrGFM_Pavg     = 43:44;
rowIbrGFM_Qavg     = 45:46;
rowIbrGFM_Angle    = 47:48;
rowIbrGFM_VRavg    = 49:52;
rowIbrGFM_VRvo     = 53:56;
rowIbrGFM_VRref    = 57:60;

% Inverters run along dim 3. The linear part has no state dependence, so the
% current guess is used only to pin down the fleet size, to check that every
% page count agrees, and to build the w0 column at the bottom of this file.
n = size(wr,3);
assert(all(cellfun(@(y)size(y,3)==n, ...
    {IphABC,vg,vo,iff,io,vi,gamma,ifRef,x,voRefDroopQ,pqavg, ...
     deltaTheta,voVRavg,voVR,voRef, ...
     Lf,rf,Cf,Rcap,Lc,rc,kc_i,kc_p,Gc,kv_i,kv_p,Gv,Mp,Mq, ...
     pset,qset,vo_q_set,wmeas_rad,w_VR_rad,Rv})))
% Every state page is a column of (re,im) slots in the layout above
assert(size(IphABC,1)==numel(idxHbIbrGFM_IphABC))
assert(size(vg,1)==numel(idxHbIbrGFM_Vg))
assert(size(vo,1)==numel(idxHbIbrGFM_Vo))
assert(size(iff,1)==numel(idxHbIbrGFM_If))
assert(size(io,1)==numel(idxHbIbrGFM_Io))
assert(size(vi,1)==numel(idxHbIbrGFM_Vi))
assert(size(gamma,1)==numel(idxHbIbrGFM_Gamma))
assert(size(ifRef,1)==numel(idxHbIbrGFM_IfRef))
assert(size(x,1)==numel(idxHbIbrGFM_X))
assert(size(wr,1)==numel(idxHbIbrGFM_Wr))
assert(size(voRefDroopQ,1)==numel(idxHbIbrGFM_VoRefDroop_q))
assert(size(pqavg,1)==numel(idxHbIbrGFM_PqAvg))
assert(size(deltaTheta,1)==numel(idxHbIbrGFM_DeltaTheta))
assert(size(voVRavg,1)==numel(idxHbIbrGFM_Vo_VR_avg))
assert(size(voVR,1)==numel(idxHbIbrGFM_Vo_VR))
assert(size(voRef,1)==numel(idxHbIbrGFM_VoRef))

% Harmonic frequency. k is an integer harmonic index, so the DC test is
% exact here.
wk = k*w0;
isDc = (k==0);
% The two restored-derivative operators, one per time base (see header).
jwkPu = 1j*wk/w0nom;        % blocks discretised against dTPU = dT*w0nom
jwkS  = 1j*wk;              % blocks discretised against dT, physical seconds
% SALT's measurement bandwidths, reproduced UNCHANGED (see header)
wmeas = wmeas_rad/(2*pi);
w_VR  = w_VR_rad/(2*pi);
eye2  = eye(2);

% Block sizes, in the order they are concatenated into JvalsLin. THIS TABLE
% IS THE INTERFACE: the caller builds its triplet map from it.
%  # | rows  | cols  | value                                      | count
%  1 |  1:6  |  1:6  | -eye(3)                 (Park I, abc side) |  36
%  2 |  7:10 |  7:10 | -eye(2)                 (Park V, dq side)  |  16
%  3 | 11:14 | 19:22 | ((-rc/Lc) - j*k*w0/w0nom)*eye2             |  16
%  4 | 11:14 | 11:14 | (1/Lc)*eye2                                |  16
%  5 | 11:14 |  7:10 | -(1/Lc)*eye2                               |  16
%  6 | 15:18 | 11:14 | -(j*k*w0/w0nom)*eye2                       |  16
%  7 | 15:18 | 15:18 | ((1/Cf) + Rcap*j*k*w0/w0nom)*eye2          |  16
%  8 | 15:18 | 19:22 | -((1/Cf) + Rcap*j*k*w0/w0nom)*eye2         |  16
%  9 | 19:22 | 15:18 | ((-rf/Lf) - j*k*w0/w0nom)*eye2             |  16
% 10 | 19:22 | 23:26 | (1/Lf)*eye2                                |  16
% 11 | 19:22 | 11:14 | -(1/Lf)*eye2                               |  16
% 12 | 23:26 | 23:26 | -eye2                                      |  16
% 13 | 23:26 | 27:30 | (kc_i + j*k*w0*kc_p)*eye2                  |  16
% 14 | 23:26 | 11:14 | Gc*eye2                                    |  16
% 15 | 27:30 | 27:30 | -(j*k*w0)*eye2                             |  16
% 16 | 27:30 | 31:34 | eye2                                       |  16
% 17 | 27:30 | 15:18 | -eye2                                      |  16
% 18 | 31:34 | 31:34 | -eye2                                      |  16
% 19 | 31:34 | 35:38 | (kv_i + j*k*w0*kv_p)*eye2                  |  16
% 20 | 31:34 | 19:22 | Gv*eye2                                    |  16
% 21 | 35:38 | 35:38 | -(j*k*w0)*eye2                             |  16
% 22 | 35:38 | 57:60 | eye2                                       |  16
% 23 | 35:38 | 11:14 | -eye2                                      |  16
% 24 | 39:40 | 39:40 | -1                                         |   4
% 25 | 39:40 | 43:44 | -Mp                                        |   4
% 26 | 41:42 | 41:42 | -1                                         |   4
% 27 | 41:42 | 45:46 | -Mq                                        |   4
% 28 | 43:44 | 43:44 | -(1 + j*k*w0/wmeas)                        |   4
% 29 | 45:46 | 45:46 | -(1 + j*k*w0/wmeas)                        |   4
% 30 | 47:48 | 47:48 | -(j*k*w0)                                  |   4
% 31 | 47:48 | 39:40 | 1                                          |   4
% 32 | 49:52 | 49:52 | (-w_VR - j*k*w0)*eye2                      |  16
% 33 | 49:52 | 19:22 | (w_VR*Rv)*eye2                             |  16
% 34 | 53:56 | 53:56 | -eye2                                      |  16
% 35 | 53:56 | 49:52 | eye2                                       |  16
% 36 | 53:56 | 19:22 | -Rv*eye2                                   |  16
% 37 | 57:60 | 57:60 | -eye2                                      |  16
% 38 | 57:60 | 53:56 | eye2                                       |  16
% 39 | 57:60 | 41:42 | [0;1]   (q row only, d row is zero)        |   8
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = [ ...
    numel(rowIbrGFM_ParkI)*numel(idxHbIbrGFM_IphABC); ...           %  1
    numel(rowIbrGFM_ParkV)*numel(idxHbIbrGFM_Vg); ...               %  2
    numel(rowIbrGFM_Coup)*numel(idxHbIbrGFM_Io); ...                %  3
    numel(rowIbrGFM_Coup)*numel(idxHbIbrGFM_Vo); ...                %  4
    numel(rowIbrGFM_Coup)*numel(idxHbIbrGFM_Vg); ...                %  5
    numel(rowIbrGFM_Cap)*numel(idxHbIbrGFM_Vo); ...                 %  6
    numel(rowIbrGFM_Cap)*numel(idxHbIbrGFM_If); ...                 %  7
    numel(rowIbrGFM_Cap)*numel(idxHbIbrGFM_Io); ...                 %  8
    numel(rowIbrGFM_Filt)*numel(idxHbIbrGFM_If); ...                %  9
    numel(rowIbrGFM_Filt)*numel(idxHbIbrGFM_Vi); ...                % 10
    numel(rowIbrGFM_Filt)*numel(idxHbIbrGFM_Vo); ...                % 11
    numel(rowIbrGFM_CurrCont)*numel(idxHbIbrGFM_Vi); ...            % 12
    numel(rowIbrGFM_CurrCont)*numel(idxHbIbrGFM_Gamma); ...         % 13
    numel(rowIbrGFM_CurrCont)*numel(idxHbIbrGFM_Vo); ...            % 14
    numel(rowIbrGFM_CurrInt)*numel(idxHbIbrGFM_Gamma); ...          % 15
    numel(rowIbrGFM_CurrInt)*numel(idxHbIbrGFM_IfRef); ...          % 16
    numel(rowIbrGFM_CurrInt)*numel(idxHbIbrGFM_If); ...             % 17
    numel(rowIbrGFM_VoltCont)*numel(idxHbIbrGFM_IfRef); ...         % 18
    numel(rowIbrGFM_VoltCont)*numel(idxHbIbrGFM_X); ...             % 19
    numel(rowIbrGFM_VoltCont)*numel(idxHbIbrGFM_Io); ...            % 20
    numel(rowIbrGFM_VoltInt)*numel(idxHbIbrGFM_X); ...              % 21
    numel(rowIbrGFM_VoltInt)*numel(idxHbIbrGFM_VoRef); ...          % 22
    numel(rowIbrGFM_VoltInt)*numel(idxHbIbrGFM_Vo); ...             % 23
    numel(rowIbrGFM_DroopW)*numel(idxHbIbrGFM_Wr); ...              % 24
    numel(rowIbrGFM_DroopW)*numel(idxHbIbrGFM_Pavg); ...            % 25
    numel(rowIbrGFM_DroopV)*numel(idxHbIbrGFM_VoRefDroop_q); ...    % 26
    numel(rowIbrGFM_DroopV)*numel(idxHbIbrGFM_Qavg); ...            % 27
    numel(rowIbrGFM_Pavg)*numel(idxHbIbrGFM_Pavg); ...              % 28
    numel(rowIbrGFM_Qavg)*numel(idxHbIbrGFM_Qavg); ...              % 29
    numel(rowIbrGFM_Angle)*numel(idxHbIbrGFM_DeltaTheta); ...       % 30
    numel(rowIbrGFM_Angle)*numel(idxHbIbrGFM_Wr); ...               % 31
    numel(rowIbrGFM_VRavg)*numel(idxHbIbrGFM_Vo_VR_avg); ...        % 32
    numel(rowIbrGFM_VRavg)*numel(idxHbIbrGFM_Io); ...               % 33
    numel(rowIbrGFM_VRvo)*numel(idxHbIbrGFM_Vo_VR); ...             % 34
    numel(rowIbrGFM_VRvo)*numel(idxHbIbrGFM_Vo_VR_avg); ...         % 35
    numel(rowIbrGFM_VRvo)*numel(idxHbIbrGFM_Io); ...                % 36
    numel(rowIbrGFM_VRref)*numel(idxHbIbrGFM_VoRef); ...            % 37
    numel(rowIbrGFM_VRref)*numel(idxHbIbrGFM_Vo_VR); ...            % 38
    numel(rowIbrGFM_VRref)*numel(idxHbIbrGFM_VoRefDroop_q); ...     % 39
    numel(rowIbrGFM_ParkI)*numel(idxHbIbrGFM_Io0)];                 % 40
assert(isequal(nBlkVals.',[36 16 16 16 16 16 16 16 16 16 16 16 16 ...
    16 16 16 16 16 16 16 16 16 16 4 4 4 4 4 4 4 4 16 16 16 16 16 16 16 8 ...
    12]))
% Running start offset of each block inside a JvalsLin column
offBlk = [0; cumsum(nBlkVals)];
nValsLin = offBlk(end);
assert(nValsLin==552)
JvalsLin = zeros(nValsLin,n);

% ------------------------- Network interface -----------------------------

% Park for I (rows 1:6)
% SALT:  Io_phA = invT*Io_dq
% HB:    -i_abc + invT2(theta)*io_dq*ibrGFMToBusI = 0
% Written abc FROM dq, so SALT's +/-120 degree balance equations disappear
% (see header (b)). Only the -i_abc half is linear; invT2(theta)*io_dq is a
% DFT sandwich and is stamped by calcHbIbrGFMMatrices_FD.
% (Block 1) d(ParkI)/d(IphABC). Constant, so the same column every inverter.
% NOTE the base change rides with the NONLINEAR half here, exactly as in
% SALT, where ibrGFMToBusI multiplies invT and IphABC enters with a bare -1.
blk_dParkI_dIphABC = -eye(numel(idxHbIbrGFM_IphABC));
JvalsLin(offBlk(1)+(1:nBlkVals(1)),:) = repmat(blk_dParkI_dIphABC(:),1,n);

% Park for V (rows 7:10)
% SALT:  V_phA = invT*Vg_dq
% HB:    -vg_dq + T2(theta)*v_abc*busToIbrGFM_V = 0
% The voltage goes the other way from the current: dq FROM abc through the
% 2x3 forward Park, which is how the TD model writes it
% ("Tv1=((T11.*va+T12.*vb)+T13.*vc).*b2V"). Only the -vg_dq half is linear.
% (Block 2) d(ParkV)/d(vg). Constant.
blk_dParkV_dVg = -eye(numel(idxHbIbrGFM_Vg));
JvalsLin(offBlk(2)+(1:nBlkVals(2)),:) = repmat(blk_dParkV_dVg(:),1,n);

% --------------------- Passive network (LCL + coupling) ------------------

% Passive coupling (rows 11:14)
% io_dq' = (1/Lc)*(-rc*io_dq + vo_dq - vg_dq) + j*w*io_dq
% => ... = 0
% HB restores the dropped io_dq' in NORMALISED time (TD: "kC=-1/dTPU+..."),
% so the operator is j*k*w0/w0nom and it enters with a MINUS, because SALT
% kept the RHS and dropped the left-hand derivative.
% j*w*io_dq is a product of two unknowns (wr and io) and goes through the
% DFT sandwich in calcHbIbrGFMMatrices_FD.
% (Block 3) d/d(io): SALT's (-rc./Lc).*eye2, plus the restored derivative.
blk_dCoup_dIo = stampReImPages(((-rc./Lc) - jwkPu).*eye2);
JvalsLin(offBlk(3)+(1:nBlkVals(3)),:) = reshape(blk_dCoup_dIo,nBlkVals(3),n);
% (Block 4) d/d(vo) and (Block 5) d/d(vg)
blk_dCoup_dVo = stampReImPages((1./Lc).*eye2);
JvalsLin(offBlk(4)+(1:nBlkVals(4)),:) = reshape(blk_dCoup_dVo,nBlkVals(4),n);
blk_dCoup_dVg = stampReImPages(-(1./Lc).*eye2);
JvalsLin(offBlk(5)+(1:nBlkVals(5)),:) = reshape(blk_dCoup_dVg,nBlkVals(5),n);

% Passive capacitor (rows 15:18)
% vo_dq' = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq + Rcap*(if_dq' - io_dq')
% =>
% 0 = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq
% TWO derivatives come back on this row, both in NORMALISED time and both
% read off the TD discretisation:
%   -(j*k*w0/w0nom)*vo_dq                from vo_dq' (TD "kP=-1/dTPU")
%   +Rcap*(j*k*w0/w0nom)*(if_dq - io_dq)  from Rcap*(if' - io')
% The second is the capacitor ESR (0.005 pu in the 39-bus data), which is
% what damps the LCL resonance; the TD model carries it as
% "gI=(1./Cf)/2+Rcp/dTPU" / "gO=(1./Cf)*(-1/2)+Rcp*(-1/dTPU)", i.e. the same
% dTPU base as the rest of the passive block and with a PLUS
% ("FP1=...+Rcp.*((ifd-ifdp)/dTPU-(iod-iodp)/dTPU)"). It falls on the if and
% io columns that already carry +/-(1/Cf), so it is ADDED INTO blocks 7 and
% 8 instead of creating two new blocks - nValsLin and the table are
% unchanged.
% (Block 6) d/d(vo): the restored vo_dq' derivative only.
blk_dCap_dVo = stampReImPages(-jwkPu.*eye2);
JvalsLin(offBlk(6)+(1:nBlkVals(6)),:) = repmat(blk_dCap_dVo(:),1,n);
% (Block 7) d/d(if) and (Block 8) d/d(io), each SALT's +/-(1/Cf) plus its
% share of the restored Rcap derivative.
capIf = (1./Cf) + Rcap.*jwkPu;
blk_dCap_dIf = stampReImPages(capIf.*eye2);
JvalsLin(offBlk(7)+(1:nBlkVals(7)),:) = reshape(blk_dCap_dIf,nBlkVals(7),n);
blk_dCap_dIo = stampReImPages(-capIf.*eye2);
JvalsLin(offBlk(8)+(1:nBlkVals(8)),:) = reshape(blk_dCap_dIo,nBlkVals(8),n);

% Passive filter (rows 19:22)
% if_dq' = (1/Lf)*(-rf*if_dq + vi_dq - vo_dq) + j*w*if_dq
% => ... = 0
% Restored in NORMALISED time (TD: "kF=-1/dTPU+..."), minus sign as above.
% (Block 9) d/d(if): SALT's (-rf./Lf).*eye2, plus the restored derivative.
blk_dFilt_dIf = stampReImPages(((-rf./Lf) - jwkPu).*eye2);
JvalsLin(offBlk(9)+(1:nBlkVals(9)),:) = reshape(blk_dFilt_dIf,nBlkVals(9),n);
% (Block 10) d/d(vi) and (Block 11) d/d(vo)
blk_dFilt_dVi = stampReImPages((1./Lf).*eye2);
JvalsLin(offBlk(10)+(1:nBlkVals(10)),:) = ...
    reshape(blk_dFilt_dVi,nBlkVals(10),n);
blk_dFilt_dVo = stampReImPages(-(1./Lf).*eye2);
JvalsLin(offBlk(11)+(1:nBlkVals(11)),:) = ...
    reshape(blk_dFilt_dVo,nBlkVals(11),n);

% --------------------------- Current controller --------------------------
% 1. vi_dq = kc_i*gamma_dq + kc_p*d(gamma_dq)/dt - j*wr*Lf*if_dq + Gc*vo_dq
% 2. d(gamma_dq)/dt = if_dq_ref - if_dq
% =>
% 1. vi_dq = kc_i*gamma_dq - j*wr*Lf*if_dq + Gc*vo_dq
% 2. 0 = if_dq_ref - if_dq
% Both derivatives are in PHYSICAL SECONDS (TD: "kG=kci/2+kcp/dT" and
% "mdt=-1/dT"), so both use j*k*w0. The first keeps the PLUS it has in TD
% ("FK1=...+kcp.*(gd-gdp)/dT+..."), the second the MINUS
% ("FKi1=-(gd-gdp)/dT+...").
% - j*wr*Lf*if_dq is a product of two unknowns -> calcHbIbrGFMMatrices_FD.
% (Block 12) d/d(vi), (Block 13) d/d(gamma) with the restored kc_p term,
% (Block 14) d/d(vo)
blk_dCurrCont_dVi = -eye(numel(idxHbIbrGFM_Vi));
JvalsLin(offBlk(12)+(1:nBlkVals(12)),:) = ...
    repmat(blk_dCurrCont_dVi(:),1,n);
blk_dCurrCont_dGamma = stampReImPages((kc_i + jwkS.*kc_p).*eye2);
JvalsLin(offBlk(13)+(1:nBlkVals(13)),:) = ...
    reshape(blk_dCurrCont_dGamma,nBlkVals(13),n);
blk_dCurrCont_dVo = stampReImPages(Gc.*eye2);
JvalsLin(offBlk(14)+(1:nBlkVals(14)),:) = ...
    reshape(blk_dCurrCont_dVo,nBlkVals(14),n);
% (Block 15) d/d(gamma) of the restored integrator derivative, (Block 16)
% d/d(if_ref), (Block 17) d/d(if)
blk_dCurrInt_dGamma = stampReImPages(-jwkS.*eye2);
JvalsLin(offBlk(15)+(1:nBlkVals(15)),:) = repmat(blk_dCurrInt_dGamma(:),1,n);
blk_dCurrInt_dIfRef = eye(numel(idxHbIbrGFM_IfRef));
JvalsLin(offBlk(16)+(1:nBlkVals(16)),:) = ...
    repmat(blk_dCurrInt_dIfRef(:),1,n);
blk_dCurrInt_dIf = -eye(numel(idxHbIbrGFM_If));
JvalsLin(offBlk(17)+(1:nBlkVals(17)),:) = repmat(blk_dCurrInt_dIf(:),1,n);

% --------------------------- Voltage controller --------------------------
% 1. if_dq_ref = kv_i*x_dq + kv_p*d(x_dq)/dt - j*wr*Cf*vo_dq + Gv*io_dq
% 2. d(x_dq)/dt = vo_dq_ref - vo_dq
% =>
% 1. if_dq_ref = kv_i*x_dq - j*wr*Cf*vo_dq + Gv*io_dq
% 2. 0 = vo_dq_ref - vo_dq
% Same time base and the same two signs as the current controller (TD:
% "kV=kvi/2+kvp/dT", "FV1c=...+kvp.*(xd-xdp)/dT+...", "FVi1=-(xd-xdp)/dT+...").
% - j*wr*Cf*vo_dq is a product of two unknowns -> calcHbIbrGFMMatrices_FD.
% (Block 18) d/d(if_ref), (Block 19) d/d(x) with the restored kv_p term,
% (Block 20) d/d(io)
blk_dVoltCont_dIfRef = -eye(numel(idxHbIbrGFM_IfRef));
JvalsLin(offBlk(18)+(1:nBlkVals(18)),:) = ...
    repmat(blk_dVoltCont_dIfRef(:),1,n);
blk_dVoltCont_dX = stampReImPages((kv_i + jwkS.*kv_p).*eye2);
JvalsLin(offBlk(19)+(1:nBlkVals(19)),:) = ...
    reshape(blk_dVoltCont_dX,nBlkVals(19),n);
blk_dVoltCont_dIo = stampReImPages(Gv.*eye2);
JvalsLin(offBlk(20)+(1:nBlkVals(20)),:) = ...
    reshape(blk_dVoltCont_dIo,nBlkVals(20),n);
% (Block 21) d/d(x) of the restored integrator derivative, (Block 22)
% d/d(vo_ref), (Block 23) d/d(vo)
blk_dVoltInt_dX = stampReImPages(-jwkS.*eye2);
JvalsLin(offBlk(21)+(1:nBlkVals(21)),:) = repmat(blk_dVoltInt_dX(:),1,n);
blk_dVoltInt_dVoRef = eye(numel(idxHbIbrGFM_VoRef));
JvalsLin(offBlk(22)+(1:nBlkVals(22)),:) = ...
    repmat(blk_dVoltInt_dVoRef(:),1,n);
blk_dVoltInt_dVo = -eye(numel(idxHbIbrGFM_Vo));
JvalsLin(offBlk(23)+(1:nBlkVals(23)),:) = repmat(blk_dVoltInt_dVo(:),1,n);

% -------------------------------- Droop ----------------------------------
% 1. wr = w0 + Mp*(pset - pavg)
% 2. vo_q_ref_droop = vo_q_set + Mq*(qset-qavg)
% 3. pavg' = wmeas*(p-pavg) s.t. p = vg_dq.' * io_dq
%    NOTE: No factor of (3/2) is needed here. It would come from the
%    Park transform for peak, line-to-neutral, non-PU quantities; with
%    RMS, line-to-line, PU quantities it cancels.
% 4. qavg' = wmeas*(q-qavg) s.t. q = vg_dq.' * [0 -1; 1 0] * io_dq
% 5. d(theta)/dt = wr
% =>
% 1. wr = w0 + Mp*(pset - pavg)
% 2. vo_q_ref_droop = vo_q_set + Mq*(qset-qavg)
% 3. 0 = p - pavg
% 4. 0 = q - qavg
% 5. ws = wr
%
% THE DROOP REFERENCE IS THE NOMINAL CONSTANT w0nom, NOT THE FREQUENCY
% UNKNOWN. SALT writes FDroopW = -wr + w0 + Mp*(pset - pavg) with w0 its
% CONSTANT property (SALT.w0 = 2*pi*SALT.f0) and keeps the solved system
% frequency in a separate variable ws, which appears only in FSpeed
% (dFspeed_dWs = -1). The TD model does the same ("FW=-wr+w0+..." with
% w0 = obj.w0). So in HB the droop reference is w0nom and the frequency
% unknown w0 enters only through the frame-angle equation below. Putting the
% unknown w0 here instead would cancel against the angle equation's wr = w0
% at DC and collapse the droop law to pavg = pset, which both removes the
% droop and leaves an islanded inverter's frequency undetermined.
% (Block 24) d/d(wr) and (Block 25) d/d(pavg)
blk_dDroopW_dWr = stampReImPages(-1);
JvalsLin(offBlk(24)+(1:nBlkVals(24)),:) = repmat(blk_dDroopW_dWr(:),1,n);
blk_dDroopW_dPavg = stampReImPages(-Mp);
JvalsLin(offBlk(25)+(1:nBlkVals(25)),:) = ...
    reshape(blk_dDroopW_dPavg,nBlkVals(25),n);
% (Block 26) d/d(vo_q_ref_droop) and (Block 27) d/d(qavg)
blk_dDroopV_dVoRefDroopQ = stampReImPages(-1);
JvalsLin(offBlk(26)+(1:nBlkVals(26)),:) = ...
    repmat(blk_dDroopV_dVoRefDroopQ(:),1,n);
blk_dDroopV_dQavg = stampReImPages(-Mq);
JvalsLin(offBlk(27)+(1:nBlkVals(27)),:) = ...
    reshape(blk_dDroopV_dQavg,nBlkVals(27),n);
% (Blocks 28,29) the two measurement filters. SALT divided these rows
% through by wmeas, so the restored derivative carries the same 1/wmeas:
% -(j*k*w0/wmeas)*pavg + p - pavg = 0, in PHYSICAL SECONDS (TD:
% "Dp1=-1/dT+wm*(-1/2)", which is this row times wmeas). p and q are
% products of two unknowns -> calcHbIbrGFMMatrices_FD.
blk_dPavg_dPavg = stampReImPages(-(1 + jwkS./wmeas));
JvalsLin(offBlk(28)+(1:nBlkVals(28)),:) = ...
    reshape(blk_dPavg_dPavg,nBlkVals(28),n);
blk_dQavg_dQavg = stampReImPages(-(1 + jwkS./wmeas));
JvalsLin(offBlk(29)+(1:nBlkVals(29)),:) = ...
    reshape(blk_dQavg_dQavg,nBlkVals(29),n);

% ------------------------------ Frame angle ------------------------------
% SALT:  FSpeed = -ws + wr, from d(theta)/dt = wr
% HB:    theta := w0*t + deltaTheta
%        => d(deltaTheta)/dt = wr - w0
%        => -j*k*w0*deltaTheta_k + wr_k - w0*delta_k0 = 0
% PHYSICAL SECONDS (TD: "Dt1=-1/dT", "Dt2=1/2"), so the operator is j*k*w0,
% and it enters with a MINUS exactly as the GENROU's block 6 does.
%
% NOTE: Should NOT replace this equation with deltaTheta=0 because
% deltaTheta at DC is just theta_DC!! Will be used accordingly in
% nonlinear equation
% (Block 30) d/d(deltaTheta). stampReIm(-1j*wk) is exactly zero at DC, so one
% branchless expression covers both cases.
blk_dAngle_dDeltaTheta = stampReImPages(-jwkS);
JvalsLin(offBlk(30)+(1:nBlkVals(30)),:) = ...
    repmat(blk_dAngle_dDeltaTheta(:),1,n);
% (Block 31) d/d(wr). Constant, stamped at every harmonic.
blk_dAngle_dWr = stampReImPages(1);
JvalsLin(offBlk(31)+(1:nBlkVals(31)),:) = repmat(blk_dAngle_dWr(:),1,n);

% --------------------------- Virtual resistor ----------------------------
% ************************
% Vo_ref == vo_q_ref_droop
% ************************
% 1. d/dt(Vo_VR_avg) = w_VR*(Rv*Io - Vo_VR_avg)
% 2. Vo_VR = Vo_VR_avg - Rv*Io
% 3. Vo_ref = vo_ref_droop + Vo_VR
% =>
% 1. 0 = w_VR*(Rv*Io - Vo_VR_avg)
% 2. Vo_VR = Vo_VR_avg - Rv*Io
% 3. Vo_ref = vo_q_ref_droop + Vo_VR
% Equation 1 is in PHYSICAL SECONDS (TD: "Ra=-1/dT-wVR/2"), and UNLIKE the
% power filters SALT kept the w_VR factor on the algebraic half, so the
% restored derivative carries NO 1/w_VR: -(j*k*w0)*vo_VR_avg + w_VR*(...).
% TD's row confirms the pairing: -1/dT against -wVR/2 and +wVR*Rv/2.
% (Block 32) d/d(vo_VR_avg) = SALT's -w_VR*eye2 plus the restored
% derivative, (Block 33) d/d(io)
blk_dVRavg_dVoVRavg = stampReImPages((-w_VR - jwkS).*eye2);
JvalsLin(offBlk(32)+(1:nBlkVals(32)),:) = ...
    reshape(blk_dVRavg_dVoVRavg,nBlkVals(32),n);
blk_dVRavg_dIo = stampReImPages((w_VR.*Rv).*eye2);
JvalsLin(offBlk(33)+(1:nBlkVals(33)),:) = ...
    reshape(blk_dVRavg_dIo,nBlkVals(33),n);
% (Blocks 34,35,36) equation 2. No derivative to restore.
blk_dVRvo_dVoVR = -eye(numel(idxHbIbrGFM_Vo_VR));
JvalsLin(offBlk(34)+(1:nBlkVals(34)),:) = repmat(blk_dVRvo_dVoVR(:),1,n);
blk_dVRvo_dVoVRavg = eye(numel(idxHbIbrGFM_Vo_VR_avg));
JvalsLin(offBlk(35)+(1:nBlkVals(35)),:) = repmat(blk_dVRvo_dVoVRavg(:),1,n);
blk_dVRvo_dIo = stampReImPages(-Rv.*eye2);
JvalsLin(offBlk(36)+(1:nBlkVals(36)),:) = ...
    reshape(blk_dVRvo_dIo,nBlkVals(36),n);
% (Blocks 37,38,39) equation 3. No derivative to restore. The droop
% reference is a q-axis quantity only, so SALT's [zeros; vo_q_ref_droop]
% becomes the 2x1 complex [0;1] - a 4x2 real block whose d rows are zero.
blk_dVRref_dVoRef = -eye(numel(idxHbIbrGFM_VoRef));
JvalsLin(offBlk(37)+(1:nBlkVals(37)),:) = repmat(blk_dVRref_dVoRef(:),1,n);
blk_dVRref_dVoVR = eye(numel(idxHbIbrGFM_Vo_VR));
JvalsLin(offBlk(38)+(1:nBlkVals(38)),:) = repmat(blk_dVRref_dVoVR(:),1,n);
blk_dVRref_dVoRefDroopQ = stampReImPages([0;1]);
JvalsLin(offBlk(39)+(1:nBlkVals(39)),:) = ...
    repmat(blk_dVRref_dVoRefDroopQ(:),1,n);

% (Block 40) d(ParkI)/d(Io0), the 0-column of the dq0 inverse Park. That
% column is ones(3,1) - no theta in it, which is exactly why it is here and
% not in the nonlinear kernel with the other two Park columns. It carries
% ibrGFMToBusI for the same reason the d and q columns do: the inverse Park
% converts device per-unit current to bus per-unit.
blk_dParkI_dIo0 = ibrGFMToBusI.*repmat(eye(2),3,1);            % [6 x 2 x n]
JvalsLin(offBlk(40)+(1:nBlkVals(40)),:) = ...
    reshape(blk_dParkI_dIo0,nBlkVals(40),n);

% ------------------------- Linear residual -------------------------------
% RESIDUAL SIGN. The caller works with J*x - F, i.e. J*x = F, so for a
% linear equation F collapses to MINUS the equation's state-free constant.
% The GENROU settles the sign: its swing constant -(KD/w)*w = -KD comes back
% as F = +KD*N. The N is there because the HB state vector holds
% UNNORMALIZED DFT coefficients (hbIdft divides by N, GAMMA does not), so a
% DC physical constant c appears as c*N - and only at k = 0, since a
% constant has no other harmonic.
%
% The GFM has exactly two state-free constants, both droop setpoints:
%   row 39 (P droop): w0nom + Mp*pset   (w0nom is the droop reference; see
%                     the Droop block above for why it is not the unknown)
%   row 41 (Q droop): vo_q_set + Mq*qset
% Both are real, so they land on the REAL row of their pair and never on the
% imaginary one. Everything else in these 60 rows has a state attached.
F = zeros(nEqns,n);
if isDc
    F(rowIbrGFM_DroopW(1),:) = -N*reshape(w0nom + Mp.*pset,1,n);
    F(rowIbrGFM_DroopV(1),:) = -N*reshape(vo_q_set + Mq.*qset,1,n);
end

% ----------------- Coupling into the system-frequency column -------------
% d/dw0 of every row at THIS harmonic. Two kinds of term contribute, and
% THEY SCALE DIFFERENTLY: w0 is a plain scalar unknown while every other
% unknown is a DFT coefficient, which at DC is N times the time-domain
% value. So a BARE ADDITIVE w0 in a DC equation needs *N - case (i) - while
% d/dw0 of something like j*k*w0*c*x does NOT, because the N is already
% inside x - case (ii).
%
% (i) The frame-angle equation's -w0, at k = 0 only. This is the GENROU's
%     dFdeltaTheta_dW0 verbatim: "wr defined by DFT, xN forces w0 to be too"
%     at DC, and "No w0 at higher harmonics by defn of fundamental" above it.
%     w0 is a real DC scalar, so it lands on the REAL row of deltaTheta.
%
% (ii) Every restored derivative, whose operator is j*k*w0 (or j*k*w0/w0nom),
%     carries a d/dw0 = j*k*(constant)*x_k onto its own row - the same
%     "every harmonic pulls on the w0 column with weight k" that the linear
%     network stamps, and the same thing stampOfHbLoadCompMatrices_FD does
%     for its stator and inertia terms. Each entry below is the derivative
%     of the corresponding block above with wk = k*w0 differentiated, times
%     that row's own harmonic coefficient. All of them vanish at k = 0, so
%     at DC only (i) survives and this reduces to SALT, where the only
%     frequency coupling is dFspeed_dWs = -1.
%
% THE COLUMN IS RETURNED IN TWO PIECES, and the split is not cosmetic - it
% is what keeps the Newton right-hand side exact. The caller adds J*x to the
% right-hand side only for NONLINEAR terms, because for a linear term J*x
% and F cancel (HB.stampGenrou). Group (i) is linear in w0, so its column
% cancels and must NOT be added; every group (ii) entry is the w0 half of a
% term BILINEAR in w0 and a state, for which J*x does not cancel and the
% contribution dFdw0*w0 has to go into the right-hand side. Hence:
%   dFdw0Bare  group (i) only   - stamp in J, keep OUT of the RHS
%   dFdw0      group (ii) only  - stamp in J, and add dFdw0*w0 to the RHS
% The Jacobian's w0 column is their sum.
dFdw0     = zeros(nEqns,n);
dFdw0Bare = zeros(nEqns,n);
% (i)
if isDc
    dFdw0Bare(rowIbrGFM_Angle(1),:) = -N;
end
% (ii) Fold the interleaved (re,im) state pages into complex coefficients
ioC         = toComplexPages(io);
voC         = toComplexPages(vo);
iffC        = toComplexPages(iff);
gammaC      = toComplexPages(gamma);
xC          = toComplexPages(x);
pqavgC      = toComplexPages(pqavg);
deltaThetaC = toComplexPages(deltaTheta);
voVRavgC    = toComplexPages(voVRavg);
% Normalised-time blocks: d/dw0 of -(j*k*w0/w0nom)*x is -(j*k/w0nom)*x
dFdw0(rowIbrGFM_Coup,:) = ...
    reshape(toReImPages(-(1j*k/w0nom).*ioC),numel(rowIbrGFM_Coup),n);
% The capacitor row restores TWO derivatives, so its w0 column carries both:
% -(j*k/w0nom)*vo from vo_dq', and +Rcap*(j*k/w0nom)*(if - io) from the ESR.
dFdw0(rowIbrGFM_Cap,:) = ...
    reshape(toReImPages((1j*k/w0nom).*(Rcap.*(iffC - ioC) - voC)), ...
    numel(rowIbrGFM_Cap),n);
dFdw0(rowIbrGFM_Filt,:) = ...
    reshape(toReImPages(-(1j*k/w0nom).*iffC),numel(rowIbrGFM_Filt),n);
% Physical-seconds blocks: d/dw0 of (+/-)(j*k*w0)*c*x is (+/-)(j*k)*c*x
dFdw0(rowIbrGFM_CurrCont,:) = ...
    reshape(toReImPages((1j*k).*kc_p.*gammaC),numel(rowIbrGFM_CurrCont),n);
dFdw0(rowIbrGFM_CurrInt,:) = ...
    reshape(toReImPages(-(1j*k).*gammaC),numel(rowIbrGFM_CurrInt),n);
dFdw0(rowIbrGFM_VoltCont,:) = ...
    reshape(toReImPages((1j*k).*kv_p.*xC),numel(rowIbrGFM_VoltCont),n);
dFdw0(rowIbrGFM_VoltInt,:) = ...
    reshape(toReImPages(-(1j*k).*xC),numel(rowIbrGFM_VoltInt),n);
dFdw0(rowIbrGFM_Pavg,:) = ...
    reshape(toReImPages(-(1j*k./wmeas).*pqavgC(1,:,:)), ...
    numel(rowIbrGFM_Pavg),n);
dFdw0(rowIbrGFM_Qavg,:) = ...
    reshape(toReImPages(-(1j*k./wmeas).*pqavgC(2,:,:)), ...
    numel(rowIbrGFM_Qavg),n);
dFdw0(rowIbrGFM_VRavg,:) = ...
    reshape(toReImPages(-(1j*k).*voVRavgC),numel(rowIbrGFM_VRavg),n);
% The frame angle takes BOTH (i) and (ii): -(j*k)*deltaTheta here, on top of
% the -N that went into dFdw0Bare. The two never coexist numerically, since
% (ii) is zero at k = 0, but they still belong to different outputs because
% only this half may reach the right-hand side.
dFdw0(rowIbrGFM_Angle,:) = ...
    reshape(toReImPages(-(1j*k).*deltaThetaC),numel(rowIbrGFM_Angle),n);
% Rows 1:10, 39:42 and 53:60 get nothing: the two Park transforms depend on
% w0 only through theta = w0*t + deltaTheta, which is pure periodicity and
% leaves the Fourier coefficients of the transformed quantity alone (the
% same argument the GENROU's nonlinear kernel makes), and the droop, VR drop
% and VR reference equations contain no w0 and no derivative.

% IphABC, vg, vi, ifRef, voVR and voRef are taken as arguments so the
% kernel's state list matches the device's variable layout one for one and
% so their page counts are checked above; no w0 coupling touches their rows.

end


function Y_ReIm = stampReImPages(Y)
%STAMPREIMPAGES  Page-aware stampReIm, i.e. [real(Y) -imag(Y); imag(Y) real(Y)]
%   written with the real and imaginary rows/cols interleaved. Identical to
%   stampReIm.m, extended with a trailing page dimension so one call covers
%   every inverter.

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
