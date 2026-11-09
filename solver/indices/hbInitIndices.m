function hbInitIndices()
%HBINITINDICES  Populate the HARMONIC-BALANCE index vectors as globals.
%
%   Call once - HB does this from its constructor - and thereafter any
%   function can reach the index vectors with a single GLOBAL declaration
%   instead of re-calling a get*Indices function on every stamp. Same
%   contract as SALTINITINDICES.
%
%   Naming: idxHb<Model>_<state>, the position of each state inside one
%   device's block AT ONE HARMONIC. The harmonic index multiplies the block
%   size; see HB.getRowIdx_fromGenrouIndices.
%
%   WHY A SEPARATE SET FROM SALTINITINDICES. These are globals, and a
%   single MATLAB session routinely holds both an HB and a SALT object
%   (run_SALT_contingency builds both), so the two index sets must not
%   collide. They also genuinely differ:
%
%     SALT    one harmonic, so the dq frame is stationary and every dq
%             quantity is a single REAL scalar. The damper currents are
%             identically zero at steady state, so SALT reduces the machine
%             to (id,iq,ifd) and drops i0,i1d,i1q,i2q.
%
%     HB      harmonics are present, so every quantity - dq included -
%             carries a (real,imag) PAIR at each harmonic, and the damper
%             currents are NOT zero: harmonic flux induces them, and the
%             dampers are the dominant harmonic path through a machine. All
%             seven windings are therefore retained.
%
%   See also SALTINITINDICES, EMTINITINDICES.

% ---------------------------------------------------------------- GENROU
% Seven stator/rotor currents, three terminal voltages, the field voltage,
% speed, rotor angle, the three-phase terminal current and the mechanical
% torque - each as a (real,imag) pair.
%
%   i = [id iq i0 | ifd i1d | i1q i2q]^T      (Kundur 3.120-3.134)
%   e = [ed eq e0]^T
global idxHbGenrou_Idq0 idxHbGenrou_Edq0 idxHbGenrou_Efd ...
       idxHbGenrou_Wr idxHbGenrou_DeltaTheta idxHbGenrou_IphABC ...
       idxHbGenrou_Tm                                           %#ok<*GVMIS>
idxHbGenrou_Idq0       =  1:14;
idxHbGenrou_Edq0       = 15:20;
idxHbGenrou_Efd        = 21:22;
idxHbGenrou_Wr         = 23:24;
% DeltaTheta, not Theta: HB writes the rotor angle as theta(t) = w0*t +
% deltaTheta(t) and carries only the deviation, because the w0*t part is
% what makes the solution periodic and is absorbed into the sample grid.
idxHbGenrou_DeltaTheta = 25:26;
idxHbGenrou_IphABC     = 27:32;
idxHbGenrou_Tm         = 33:34;

% --------------------------------------------------------------- EXCDC1
% 5 controller states + 1 terminal-voltage calculation, as (real,imag).
global idxHbExcDc1_x1 idxHbExcDc1_x2 idxHbExcDc1_x3 idxHbExcDc1_x4 ...
       idxHbExcDc1_x5 idxHbExcDc1_Vt
idxHbExcDc1_x1 =  1:2;   idxHbExcDc1_x2 =  3:4;   idxHbExcDc1_x3 =  5:6;
idxHbExcDc1_x4 =  7:8;   idxHbExcDc1_x5 =  9:10;  idxHbExcDc1_Vt = 11:12;

% --------------------------------------------------------------- TGOV1D
% 2 states + 1 mechanical-power calculation, as (real,imag).
global idxHbTgov1d_x1 idxHbTgov1d_x2 idxHbTgov1d_Pmech
idxHbTgov1d_x1 = 1:2;  idxHbTgov1d_x2 = 3:4;  idxHbTgov1d_Pmech = 5:6;

% ------------------------------------------- COMPOSITE LOAD (IM + ZIP)
% SALT carries 20 scalars here because its dq quantities are real DC
% values; HB needs a (real,imag) pair for each at every harmonic, so a dq
% pair costs 4 slots and a scalar 2. The abc terminal current was already a
% phasor triple in SALT and stays 6.
%
% BALANCED, MEANING ZERO-SUM - NOT positive-sequence-only. SALT ties the
% three phases together with iphBC = phShiftBC*iphA, i.e. i_b = a*i_a and
% i_c = a^2*i_a. That is the POSITIVE-SEQUENCE condition: two complex
% constraints. It implies zero sum (1 + a + a^2 = 0) and is strictly
% stronger, and at a single fundamental harmonic with a balanced source and
% network the two coincide - which is why SALT is right as written.
%
% With more than one harmonic they part company, and the condition the
% model actually rests on is the weaker one. Lecture 7 of the circuit
% simulator course states it directly for the TPIM: "assume balanced system
% and ignore the zero frame", which in the time domain is
%
%     i_a + i_b + i_c = 0
%
% - ONE scalar constraint, one complex equation per harmonic. HB uses that.
% The counts then close exactly (34 equations for 34 unknowns here, 60 for
% 60 in the IBR) with no zero-sequence inductance invented anywhere,
% because the zero frame is DROPPED rather than modelled.
%
% What that excludes is zero-sequence, i.e. triplen harmonics (3rd, 9th,
% ...). Negative sequence is fine - zero-sum admits it - so a 2nd harmonic,
% which is normally negative sequence, is representable. Lifting SALT's
% +/-120 degree form naively would have silently forced it positive.
%
% In practice the constraint never appears as its own equation: the
% terminal current is written as abc FROM dq through the 3x2 inverse Park,
% and every column of that matrix is a cos or a sin evaluated at the three
% phase angles, which already sums to zero - so zero sum holds by
% construction. The terminal voltage goes the other way, dq FROM abc
% through the 2x3 forward Park. That split mirrors the lecture's own
% circuit: dependent voltage sources feeding the machine, dependent current
% sources feeding the network.
%
% The SIGN of the sin in those matrices is not the same for the two
% devices, and it is not a free choice either way: it is pinned by the
% speed-voltage matrix each model already carries. The composite load uses
% the -sin convention, so its dq pair IS the phasor and SALT's
% wMat = [0 -w; w 0] carries over verbatim; the GFM inverter uses +sin
% (getParkMatrices_2x3_PSCAD), so its dq pair is the phasor's CONJUGATE and
% SALT's passive equations carry -stampReIm(jw) instead. Each kernel's
% header derives its own; zero sum holds in both.
global idxHbLoadComp_Vs idxHbLoadComp_IphABC idxHbLoadComp_WSlip ...
       idxHbLoadComp_DeltaThetaSys idxHbLoadComp_Is idxHbLoadComp_Ir ...
       idxHbLoadComp_Wr idxHbLoadComp_PQ idxHbLoadComp_Vm ...
       idxHbLoadComp_Izip
% Induction motor
idxHbLoadComp_Vs            =  1:4;    % (vd_s,vq_s)
idxHbLoadComp_IphABC        =  5:10;   % three phases
idxHbLoadComp_WSlip         = 11:12;
% DeltaThetaSys, not ThetaSys: the system frame angle is theta_sys =
% w0*t + deltaThetaSys, and only the deviation is carried, for the same
% reason the GENROU carries deltaTheta.
idxHbLoadComp_DeltaThetaSys = 13:14;
idxHbLoadComp_Is            = 15:18;   % (id_s,iq_s)
idxHbLoadComp_Ir            = 19:22;   % (id_r,iq_r)
idxHbLoadComp_Wr            = 23:24;
% ZIP
idxHbLoadComp_PQ            = 25:28;   % (P,Q)
idxHbLoadComp_Vm            = 29:30;
idxHbLoadComp_Izip          = 31:34;   % (id,iq)

% ------------------------------------------------------- GFM IBR (grid
% forming). SALT's 33 scalars become 60 for the same reason as above.
global idxHbIbrGFM_IphABC idxHbIbrGFM_Vg idxHbIbrGFM_Vo idxHbIbrGFM_If ...
       idxHbIbrGFM_Io idxHbIbrGFM_Vi idxHbIbrGFM_Gamma ...
       idxHbIbrGFM_IfRef idxHbIbrGFM_X idxHbIbrGFM_Wr ...
       idxHbIbrGFM_VoRefDroop_q idxHbIbrGFM_PqAvg ...
       idxHbIbrGFM_DeltaTheta idxHbIbrGFM_Vo_VR_avg idxHbIbrGFM_Vo_VR ...
       idxHbIbrGFM_VoRef idxHbIbrGFM_Io0
% Network interface
idxHbIbrGFM_IphABC       =  1:6;
% Passive network (LCL filter + coupling)
idxHbIbrGFM_Vg           =  7:10;
idxHbIbrGFM_Vo           = 11:14;
idxHbIbrGFM_If           = 15:18;
idxHbIbrGFM_Io           = 19:22;
idxHbIbrGFM_Vi           = 23:26;
% Current controller
idxHbIbrGFM_Gamma        = 27:30;
% Voltage controller
idxHbIbrGFM_IfRef        = 31:34;
idxHbIbrGFM_X            = 35:38;
% Droop
idxHbIbrGFM_Wr           = 39:40;
idxHbIbrGFM_VoRefDroop_q = 41:42;
idxHbIbrGFM_PqAvg        = 43:46;
% DeltaTheta, not Theta: same convention as the GENROU and the composite
% load. SALT's "d(theta)/dt = wr" becomes "d(deltaTheta)/dt = wr - w0",
% which at DC pins wr to w0 and at higher harmonics reads
% j*k*w0*deltaTheta_k = wr_k.
idxHbIbrGFM_DeltaTheta   = 47:48;
% Virtual resistor
idxHbIbrGFM_Vo_VR_avg    = 49:52;
idxHbIbrGFM_Vo_VR        = 53:56;
idxHbIbrGFM_VoRef        = 57:60;
% Zero-sequence current of the dq0 inverse Park. Appended rather than
% inserted so none of the indices above move. See the DELTA, NOT WYE note
% at the top of CALCHBIBRGFMMATRICES_FD.
idxHbIbrGFM_Io0          = 61:62;

end
