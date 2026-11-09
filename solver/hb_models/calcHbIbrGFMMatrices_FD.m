function [JvalsNlin,Fnlin,JvalsVabc] = calcHbIbrGFMMatrices_FD( ...
    ... % Time-domain previous guess (pages: [nVar x N x numIbrGFM])
    vgTd,voTd,ifTd,ioTd,wrTd,deltaThetaTd, vabcTd, ...
    ... % Frequency-domain previous guess, used for J*x
    w0, ...
    ... % Parameters (pages: [1 x 1 x numIbrGFM])
    Lf,Cf, w0nom, ...
    ... % Transform
    ibrGFMToBusI,busToIbrGFM_V, K,N,GAMMA,invGAMMA,tSamp)
%CALCHBIBRGFMMATRICES_FD  Nonlinear harmonic-balance GFM-IBR kernel.
%
%   DELTA, NOT WYE - WHY THIS MODEL IS NOT IDENTICAL TO SALT.
%
%   This GFM-IBR model is not identical to SALT, since it puts the
%   zero-sequence term on voltage (i.e. va+vb+vc=0). This is because the
%   GFM-IBR is a voltage source. The equivalent circuit requires a delta
%   connection on the network side (i.e. va+vb+vc=0) to set the voltages
%   instead of a wye connection (i.e. ia+ib+ic=0) to set the currents;
%   otherwise, we will have two voltage sources connected to each other in
%   series (i.e. voltage source of IBR connected with controlled voltage
%   source representing dq-transform) which is illegal. This is not a
%   problem in SALT because there is not a "DC" harmonic where all
%   capacitor shunts are open, but in HB, the DC harmonic makes the source
%   floating on the IBR device side, so the two voltage sources are truly
%   open.
%
%   Note that despite these differences between HB and SALT, the results
%   should still be the same because SALT assumes balanced, which means
%   enforcing va+vb+vc=0 is equivalent to ia+ib+ic=0.
%
%   WHERE IT LIVES. The 2x3 Park pair became a genuine 3x3 dq0 pair, which
%   costs the device two rows (60 -> 62) and splits across three places,
%   none of them here, because neither new piece depends on theta:
%     inverse Park 0-column, ones(3,1)   the io0 the inverter draws, a
%                                        constant times ibrGFMToBusI
%                                        -> block 40 of the LINEAR kernel
%     forward Park 0-row, (1/3)(1,1,1)   the equation v0 = 0 against the
%                                        terminal bus nodes
%                                        -> the ibrGFMZeroSeq slice in HB
%     io0 variable, ParkZero equation    rows and columns 61:62
%                                        -> HBINITINDICES, numIbrGFMEqns
%   Under balance io0 solves to zero, which is the numerical form of the
%   last paragraph above and is asserted in validation/validateHbVsSalt.
%
%   Dense-across-harmonics part of the grid-forming inverter equations: each
%   derivative is formed in the time domain and squeezed back into harmonic
%   coefficients by the DFT sandwich (GAMMA .* dRow) * invGAMMA. Every input
%   carries the inverters along the third dimension, so the whole fleet is
%   evaluated in one pass. JvalsNlin is [nValsNlin x numIbrGFM] with the
%   dense blocks flattened column-major and concatenated in the fixed order
%   documented below; Fnlin is the nonlinear residual [62*(K+1) x numIbrGFM];
%   JvalsVabc is the one block that leaves the device, coupling the Park V
%   rows to the terminal bus voltage.
%
%   Generalised from stampOfIbrGFMMatrices_SS.m. The linear, harmonic-
%   diagonal half - including every restored d/dt - is in
%   stampOfHbIbrGFMMatrices_FD; read its header first, because it carries the
%   time-base audit, the Park sign convention, the droop reference and the
%   two derivatives restored on the capacitor row (vo' and the Rcap ESR).
%   Nothing here needs Rcap: that term is a constant times (if - io), hence
%   linear and harmonic-diagonal.
%
%   WHAT IS NONLINEAR HERE, AND WHY EACH ONE IS:
%
%     Park I, Park V   trig of an UNKNOWN angle times an UNKNOWN. SALT
%                      short-circuits both with a 2x2 invT onto phase A,
%                      which only holds at one harmonic; HB uses the genuine
%                      3x2 inverse Park for the current (dq -> abc) and the
%                      2x3 forward Park for the voltage (abc -> dq), with
%                      theta(t) = w0*t + deltaTheta(t).
%
%     every j*w*x      SALT writes jw = stamp_j_invW0.*wr, i.e. the PRODUCT
%                      of the harmonic-expanded frame speed wr with another
%                      unknown. Five of them: j*w*io (coupling), j*w*vo
%                      (capacitor), j*w*if (filter), -j*wr*Lf*if (current
%                      controller) and -j*wr*Cf*vo (voltage controller).
%                      The 1/w0nom inside stamp_j_invW0 is the per-unit
%                      inductance's normalisation, NOT a frequency: it is the
%                      same constant the TD model spells "jw=1/w0" and uses
%                      with dTPU, and it does NOT become the unknown w0.
%
%     p and q          p = vg_dq.'*io_dq and q = vg_dq.'*[0 -1;1 0]*io_dq,
%                      bilinear in two unknowns.
%
%   PARK SIGN CONVENTION: getParkMatrices_2x3_PSCAD, which the TD model
%   inlines, gives T = 2/3*[cos(theta_abc); sin(theta_abc)] and
%   JT = 2/3*[-sin(theta_abc); cos(theta_abc)] with
%   theta_abc = theta + [0 -2*pi/3 2*pi/3]. Its matching 3x2 inverse is
%   invT2 = [cos(theta_abc).', +sin(theta_abc).'], which satisfies
%   T*invT2 = eye(2). That +sin DIFFERS from the GENROU's
%   getInvParkMatrices ([cos, -sin, 1]); this device keeps the PSCAD
%   convention, as SALT's own 2x2 invT = [cos sin; sin -cos] already did.
%   Because every column of invT2 sums to zero over the three phases,
%   i_a+i_b+i_c = 0 holds by construction and SALT's FIphBC disappears.
%
%   See also STAMPOFHBIBRGFMMATRICES_FD, CALCHBGENROUMATRICES_FD,
%   GETPARKMATRICES_2X3_PSCAD, STAMPOFIBRGFMMATRICES_SS.

% VARIABLE LAYOUT (fixed): one GFM-IBR block per harmonic is 62 rows, every
% quantity a (real,imag) pair, so a global row is k*62 + localIndex, i.e.
% harmonic-outer / variable-inner, and a dq pair is stored
% (d_re, d_im, q_re, q_im). That is exactly the ordering both the DFT
% sandwich and hbDft produce, so no block needs reordering.
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
% Io0                61:62      dq0 zero-sequence current (see header)
idxHbIbrGFM_Vg         =  7:10;
idxHbIbrGFM_Vo         = 11:14;
idxHbIbrGFM_If         = 15:18;
idxHbIbrGFM_Io         = 19:22;
idxHbIbrGFM_Wr         = 39:40;
idxHbIbrGFM_DeltaTheta = 47:48;
nEqns = 62;
% Equation rows touched by this kernel (the full row table is in
% stampOfHbIbrGFMMatrices_FD; note it is NOT the same permutation as the
% variables, e.g. the io equation is on rows 11:14 while io is on 19:22).
rowIbrGFM_ParkI    =  1:6;
rowIbrGFM_ParkV    =  7:10;
rowIbrGFM_Coup     = 11:14;
rowIbrGFM_Cap      = 15:18;
rowIbrGFM_Filt     = 19:22;
rowIbrGFM_CurrCont = 23:26;
rowIbrGFM_VoltCont = 31:34;
rowIbrGFM_Pavg     = 43:44;
rowIbrGFM_Qavg     = 45:46;
% The terminal bus voltage is three nodes, (re,im) each, at the SAME harmonic
nVarsAbc = 3;

% K is # harmonics, NOT including zeros, so total is +1
K_tot = K + 1;
% Inverters run along dim 3
n = size(wrTd,3);
% Time runs along dim 2; the scalar states carry one row and tSamp is a
% column of samples (same contract as the GENROU kernel).
assert(size(wrTd,1)==1 && size(deltaThetaTd,1)==1 && size(tSamp,2)==1)
assert(size(vgTd,1)==2 && size(voTd,1)==2 && size(ifTd,1)==2 && ...
    size(ioTd,1)==2 && size(vabcTd,1)==nVarsAbc)
assert(size(vgTd,2)==N && size(voTd,2)==N && size(ifTd,2)==N && ...
    size(ioTd,2)==N && size(vabcTd,2)==N && size(wrTd,2)==N && ...
    size(deltaThetaTd,2)==N && size(tSamp,1)==N)
assert(all(cellfun(@(y)size(y,3)==n, ...
    {vgTd,voTd,ifTd,ioTd,deltaThetaTd,vabcTd, ...
     Lf,Cf,ibrGFMToBusI,busToIbrGFM_V})))
% Each DFT sandwich must come out [2*K_tot x 2*K_tot] to land on a harmonic
% block, which holds only when GAMMA keeps exactly DC..K (N = 2*K+1).
assert(size(GAMMA,1)==2*K_tot && size(GAMMA,2)==N)
assert(size(invGAMMA,1)==N && size(invGAMMA,2)==2*K_tot)

% Block layout, in the order they are concatenated into JvalsNlin. Each is a
% dense (numel(rows)*K_tot) x (numel(cols)*K_tot) array, rows
% reshape((0:K)*62 + rows(:),[],1) and cols likewise. THIS TABLE IS THE
% INTERFACE: the caller builds its triplet map from it.
%  # | rows  | cols  | term                                  | count/K_tot^2
%  1 |  1:6  | 19:22 | d(invT2(theta)*io)/d(io)              |  24
%  2 |  1:6  | 47:48 | d(invT2(theta)*io)/d(deltaTheta)      |  12
%  3 |  7:10 | 47:48 | d(T2(theta)*v_abc)/d(deltaTheta)      |   8
%  4 | 11:14 | 19:22 | d(j*w*io)/d(io)                       |  16
%  5 | 11:14 | 39:40 | d(j*w*io)/d(wr)                       |   8
%  6 | 15:18 | 11:14 | d(j*w*vo)/d(vo)                       |  16
%  7 | 15:18 | 39:40 | d(j*w*vo)/d(wr)                       |   8
%  8 | 19:22 | 15:18 | d(j*w*if)/d(if)                       |  16
%  9 | 19:22 | 39:40 | d(j*w*if)/d(wr)                       |   8
% 10 | 23:26 | 15:18 | d(-j*wr*Lf*if)/d(if)                  |  16
% 11 | 23:26 | 39:40 | d(-j*wr*Lf*if)/d(wr)                  |   8
% 12 | 31:34 | 11:14 | d(-j*wr*Cf*vo)/d(vo)                  |  16
% 13 | 31:34 | 39:40 | d(-j*wr*Cf*vo)/d(wr)                  |   8
% 14 | 43:44 |  7:10 | d(p)/d(vg)                            |   8
% 15 | 43:44 | 19:22 | d(p)/d(io)                            |   8
% 16 | 45:46 |  7:10 | d(q)/d(vg)                            |   8
% 17 | 45:46 | 19:22 | d(q)/d(io)                            |   8
% Returned separately, NOT part of JvalsNlin:
%    |  7:10 | bus   | d(T2(theta)*v_abc)/d(v_abc)           |  24
% The two dq0 zero-sequence blocks are NOT here and never will be: both are
% theta-free, so rows 1:6 x cols 61:62 is block 40 of the linear kernel and
% rows 61:62 x bus is the constant ibrGFMZeroSeq slice. See the header.
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = K_tot^2 * [ ...
    numel(rowIbrGFM_ParkI)*numel(idxHbIbrGFM_Io); ...               %  1
    numel(rowIbrGFM_ParkI)*numel(idxHbIbrGFM_DeltaTheta); ...       %  2
    numel(rowIbrGFM_ParkV)*numel(idxHbIbrGFM_DeltaTheta); ...       %  3
    numel(rowIbrGFM_Coup)*numel(idxHbIbrGFM_Io); ...                %  4
    numel(rowIbrGFM_Coup)*numel(idxHbIbrGFM_Wr); ...                %  5
    numel(rowIbrGFM_Cap)*numel(idxHbIbrGFM_Vo); ...                 %  6
    numel(rowIbrGFM_Cap)*numel(idxHbIbrGFM_Wr); ...                 %  7
    numel(rowIbrGFM_Filt)*numel(idxHbIbrGFM_If); ...                %  8
    numel(rowIbrGFM_Filt)*numel(idxHbIbrGFM_Wr); ...                %  9
    numel(rowIbrGFM_CurrCont)*numel(idxHbIbrGFM_If); ...            % 10
    numel(rowIbrGFM_CurrCont)*numel(idxHbIbrGFM_Wr); ...            % 11
    numel(rowIbrGFM_VoltCont)*numel(idxHbIbrGFM_Vo); ...            % 12
    numel(rowIbrGFM_VoltCont)*numel(idxHbIbrGFM_Wr); ...            % 13
    numel(rowIbrGFM_Pavg)*numel(idxHbIbrGFM_Vg); ...                % 14
    numel(rowIbrGFM_Pavg)*numel(idxHbIbrGFM_Io); ...                % 15
    numel(rowIbrGFM_Qavg)*numel(idxHbIbrGFM_Vg); ...                % 16
    numel(rowIbrGFM_Qavg)*numel(idxHbIbrGFM_Io)];                   % 17
assert(isequal(nBlkVals.',K_tot^2*[24 12 8 16 8 16 8 16 8 16 8 16 8 ...
    8 8 8 8]))
% Running start offset of each block inside a JvalsNlin column
offBlk = [0; cumsum(nBlkVals)];
nValsNlin = offBlk(end);
assert(nValsNlin==196*K_tot^2)
JvalsNlin = zeros(nValsNlin,n);

% ------------------------------ Park angle -------------------------------
% NOTE: Need to double-check!! Since w0 is fundamental frequency, there
% should be a derivative term w.r.t. w0. However, its role in the Park
% transform is purely enforcing periodicty (e.g. cos(w0*t + deltaTheta)) so
% in the frequency domain, the change of w0 doesn't affect the fourier
% coefficients of the outcome of the Park transform...
%
% tSamp is transposed so time lines up with dim 2 of the page arrays
thetaTd = w0*(tSamp.') + deltaThetaTd;              % [1 x N x n]
% Inlined from getParkMatrices_2x3_PSCAD so the phase offsets can carry a
% page dimension: that helper builds a 2x3 for a row vector of angles, which
% cannot also hold pages.
thetaTd_phAbc = thetaTd + [0;-2*pi/3;2*pi/3];       % [3 x N x n]
cosTd = cos(thetaTd_phAbc);
sinTd = sin(thetaTd_phAbc);

% Split the dq pages once; (1,:,:) is d and (2,:,:) is q
vgdTd = vgTd(1,:,:);  vgqTd = vgTd(2,:,:);
vodTd = voTd(1,:,:);  voqTd = voTd(2,:,:);
ifdTd = ifTd(1,:,:);  ifqTd = ifTd(2,:,:);
iodTd = ioTd(1,:,:);  ioqTd = ioTd(2,:,:);

% ---------------------- Park for I (rows 1:6) ----------------------------
% SALT:  Io_phA = invT*Io_dq
% HB:    -i_abc + invT2(theta)*io_dq*ibrGFMToBusI = 0
% invT2 row p = [cos(theta_p)  sin(theta_p)]   (PSCAD +sin, see header)
% d(invT2)/dtheta row p = [-sin(theta_p)  cos(theta_p)]
% ibrGFMToBusI rides with this half exactly as it does in SALT, where it
% multiplies invT and IphABC enters the linear stamp with a bare -1.
% (Val)
InvParkI = ibrGFMToBusI.*(cosTd.*iodTd + sinTd.*ioqTd);         % [3 x N x n]
% (Deriv)
dInvParkI_dIod        = ibrGFMToBusI.*cosTd;
dInvParkI_dIoq        = ibrGFMToBusI.*sinTd;
dInvParkI_dDeltaTheta = ibrGFMToBusI.*(cosTd.*ioqTd - sinTd.*iodTd);

% ---------------------- Park for V (rows 7:10) ---------------------------
% SALT:  V_phA = invT*Vg_dq
% HB:    -vg_dq + T2(theta)*v_abc*busToIbrGFM_V = 0
% The voltage goes the opposite way from the current: dq FROM abc through the
% 2x3 forward Park, which is how the TD model writes it
% ("Tv1=((T11.*va+T12.*vb)+T13.*vc).*b2V" with T11=2/3*cos(theta_a)).
% Rows of Tcos/Tsin are the three phases, i.e. the two rows of T2 laid out
% per phase, with the base change folded in.
Tcos = (2/3)*busToIbrGFM_V.*cosTd;                  % [3 x N x n]
Tsin = (2/3)*busToIbrGFM_V.*sinTd;
% (Val) contract over the three phases at each time sample
ParkV = [sum(Tcos.*vabcTd,1); sum(Tsin.*vabcTd,1)];             % [2 x N x n]
% (Deriv) d(T2)/dtheta = 2/3*[-sin; cos], i.e. JT of getParkMatrices_2x3_PSCAD
dParkV_dDeltaTheta = [sum(-Tsin.*vabcTd,1); sum(Tcos.*vabcTd,1)];
% ...and w.r.t. the bus voltage itself, which is the block that leaves the
% device. Tcos/Tsin ARE those derivatives, one per phase.

% --------------- Speed-voltage products, j*w*x (rows 11:34) --------------
% jw = stamp_j_invW0.*wr = (1/w0nom)*[0 -wr; wr 0], so
%   jw*x = (1/w0nom)*[-wr*xq; wr*xd]
% and d(jw*x)/dx is off-diagonal only, d(jw*x)/dwr is (1/w0nom)*[-xq; xd].
% The 1/w0nom is the per-unit inductance normalisation that goes with the
% dTPU time base, NOT a frequency - see the header.
% (Val)
jwIo = [-wrTd.*ioqTd; wrTd.*iodTd]./w0nom;          % coupling, rows 11:14
jwVo = [-wrTd.*voqTd; wrTd.*vodTd]./w0nom;          % capacitor, rows 15:18
jwIf = [-wrTd.*ifqTd; wrTd.*ifdTd]./w0nom;          % filter, rows 19:22
% (Deriv) d(jw*x)/dx is the same row for all three, so the sandwich is
% formed once and the block reused - as the GENROU reuses dInvParkV_dV for
% the current.
dJw_dXq =  -wrTd./w0nom;                            % d(d-row)/d(x_q)
dJw_dXd =   wrTd./w0nom;                            % d(q-row)/d(x_d)

% Current controller (rows 23:26)
% vi_dq = kc_i*gamma_dq + kc_p*d(gamma_dq)/dt - j*wr*Lf*if_dq + Gc*vo_dq
% The nonlinear term is -Lf*(jw*if) = (Lf/w0nom)*[wr*ifq; -wr*ifd], which is
% the TD model's "+(mL.*sK2)/2" / "+(-mL.*sK1)/2" with mL = Lf/w0nom.
LfPu = Lf./w0nom;
currContJw = [LfPu.*(wrTd.*ifqTd); -LfPu.*(wrTd.*ifdTd)];
dCurrContJw_dIfq =  LfPu.*wrTd;
dCurrContJw_dIfd = -LfPu.*wrTd;
dCurrContJw_dWr  = [LfPu.*ifqTd; -LfPu.*ifdTd];

% Voltage controller (rows 31:34)
% if_dq_ref = kv_i*x_dq + kv_p*d(x_dq)/dt - j*wr*Cf*vo_dq + Gv*io_dq
% Nonlinear term -Cf*(jw*vo) = (Cf/w0nom)*[wr*voq; -wr*vod], the TD model's
% "+(mC.*sV2)/2" / "+(-mC.*sV1)/2" with mC = Cf/w0nom.
CfPu = Cf./w0nom;
voltContJw = [CfPu.*(wrTd.*voqTd); -CfPu.*(wrTd.*vodTd)];
dVoltContJw_dVoq =  CfPu.*wrTd;
dVoltContJw_dVod = -CfPu.*wrTd;
dVoltContJw_dWr  = [CfPu.*voqTd; -CfPu.*vodTd];

% ------------------- Measured power (rows 43:46) -------------------------
% p = vg_dq.' * io_dq
% q = vg_dq.' * [0 -1; 1 0] * io_dq
%    NOTE: No factor of (3/2) is needed here. It would come from the
%    Park transform for peak, line-to-neutral, non-PU quantities; with
%    RMS, line-to-line, PU quantities it cancels.
% (Val) contract dim 1 per time sample rather than calling pagemtimes with a
% time-indexed middle dimension
pTd = vgdTd.*iodTd + vgqTd.*ioqTd;
qTd = vgqTd.*iodTd - vgdTd.*ioqTd;      % TD: "vgjio=vgd.*(-ioq)+vgq.*iod"
% (Deriv) both are bilinear, so both factors get a block, exactly as the
% GENROU's quadratic torque does

% ------------------------- Jacobian values -------------------------------
% The 17 blocks are built in the documented order and flattened in one pass,
% so the offsets cannot drift from the table above. sand() is NLINBLOCK with
% the transform arguments bound, purely to keep these lines readable.
sand = @(dRows) nlinBlock(GAMMA,invGAMMA,K,n,dRows);
% d(jw*x)/dx is the SAME derivative row for the coupling, the capacitor and
% the filter, so the sandwich is formed once and the block reused - as the
% GENROU reuses dInvParkV_dV for its current Park.
blk_dJw_dX = sand({[],dJw_dXq; dJw_dXd,[]});
blk = { ...
    sand([rowCell(dInvParkI_dIod), rowCell(dInvParkI_dIoq)]); ...   %  1
    sand(rowCell(dInvParkI_dDeltaTheta)); ...                       %  2
    sand(rowCell(dParkV_dDeltaTheta)); ...                          %  3
    blk_dJw_dX; ...                                                 %  4
    sand({-ioqTd./w0nom; iodTd./w0nom}); ...                        %  5
    blk_dJw_dX; ...                                                 %  6
    sand({-voqTd./w0nom; vodTd./w0nom}); ...                        %  7
    blk_dJw_dX; ...                                                 %  8
    sand({-ifqTd./w0nom; ifdTd./w0nom}); ...                        %  9
    sand({[],dCurrContJw_dIfq; dCurrContJw_dIfd,[]}); ...           % 10
    sand(rowCell(dCurrContJw_dWr)); ...                             % 11
    sand({[],dVoltContJw_dVoq; dVoltContJw_dVod,[]}); ...           % 12
    sand(rowCell(dVoltContJw_dWr)); ...                             % 13
    sand({iodTd, ioqTd}); ...                                       % 14
    sand({vgdTd, vgqTd}); ...                                       % 15
    sand({-ioqTd, iodTd}); ...                                      % 16
    sand({vgqTd, -vgdTd})};                                         % 17
assert(numel(blk)==numel(nBlkVals))
for idxBlk = 1:numel(blk)
    JvalsNlin(offBlk(idxBlk)+(1:nBlkVals(idxBlk)),:) = ...
        reshape(blk{idxBlk},nBlkVals(idxBlk),n);
end

% --------------- Park V coupling into the terminal bus voltage -----------
% Rows 7:10 against the three bus nodes at every harmonic pair, ordered
% harmonic-outer / node-inner, i.e. [4*K_tot x 6*K_tot]. For the GENROU this
% coupling is an identity and therefore linear; here it is T2(theta), so the
% whole thing is a DFT sandwich and belongs in this kernel.
JvalsVabc = sand([rowCell(Tcos).'; rowCell(Tsin).']);
assert(size(JvalsVabc,1)==numel(rowIbrGFM_ParkV)*K_tot)
assert(size(JvalsVabc,2)==2*nVarsAbc*K_tot)
assert(size(JvalsVabc,3)==n)
JvalsVabc = reshape(JvalsVabc, ...
    numel(rowIbrGFM_ParkV)*K_tot*2*nVarsAbc*K_tot, n);

% ----------------------- Nonlinear residual ------------------------------
% Fnlin holds the harmonic coefficients of the nonlinear terms themselves,
% g(x_prev), at the rows the terms belong to, with the PLAIN sign of the
% term. The caller forms F = JvalsNlin*x_prev - Fnlin (GENROU convention), so
% the minus on g(x_prev) lives in the caller.
Fnlin = zeros(nEqns*K_tot,n);
Fnlin(hbRows(rowIbrGFM_ParkI,K,nEqns),:)    = hbDftPages(GAMMA,K,InvParkI);
Fnlin(hbRows(rowIbrGFM_ParkV,K,nEqns),:)    = hbDftPages(GAMMA,K,ParkV);
Fnlin(hbRows(rowIbrGFM_Coup,K,nEqns),:)     = hbDftPages(GAMMA,K,jwIo);
Fnlin(hbRows(rowIbrGFM_Cap,K,nEqns),:)      = hbDftPages(GAMMA,K,jwVo);
Fnlin(hbRows(rowIbrGFM_Filt,K,nEqns),:)     = hbDftPages(GAMMA,K,jwIf);
Fnlin(hbRows(rowIbrGFM_CurrCont,K,nEqns),:) = hbDftPages(GAMMA,K,currContJw);
Fnlin(hbRows(rowIbrGFM_VoltCont,K,nEqns),:) = hbDftPages(GAMMA,K,voltContJw);
Fnlin(hbRows(rowIbrGFM_Pavg,K,nEqns),:)     = hbDftPages(GAMMA,K,pTd);
Fnlin(hbRows(rowIbrGFM_Qavg,K,nEqns),:)     = hbDftPages(GAMMA,K,qTd);

end


function B = nlinBlock(GAMMA,invGAMMA,K,n,dRows)
%NLINBLOCK  One dense harmonic block from a table of time-domain derivatives.
%   dRows is an [nRowVar x nColVar] cell of [1 x N x nDev] derivative rows,
%   one per (equation, variable) complex pair; an EMPTY entry is a structural
%   zero and is skipped rather than sandwiched. The result is
%   [2*nRowVar*(K+1) x 2*nColVar*(K+1) x nDev], rows and cols both ordered
%   harmonic-outer / variable-inner with (re,im) innermost, which is the
%   block ordering the caller's triplet map expects.

[nRowVar,nColVar] = size(dRows);
K_tot = K + 1;
B = zeros(2*nRowVar*K_tot,2*nColVar*K_tot,n);

for idxVar = 1:nRowVar
    for jdxVar = 1:nColVar
        if isempty(dRows{idxVar,jdxVar})
            continue
        end
        assert(size(dRows{idxVar,jdxVar},3)==n)
        B(idxTermRow_ReIm(idxVar,nRowVar,K), ...
            idxTermRow_ReIm(jdxVar,nColVar,K),:) = ...
            dftSandwich(GAMMA,invGAMMA,dRows{idxVar,jdxVar});
    end
end

end


function C = rowCell(A)
%ROWCELL  Split a [nRow x N x nDev] derivative array into an nRow x 1 cell of
%   [1 x N x nDev] rows, i.e. the single-column form NLINBLOCK takes.

C = num2cell(A,[2 3]);

end


function S = dftSandwich(GAMMA,invGAMMA,dRow)
%DFTSANDWICH  (GAMMA .* dRow) * invGAMMA for every inverter at once.
%   dRow is [1 x N x nDev]: one time-domain derivative row per inverter.
%   The broadcast ".*" IS the diagonal multiplication by the derivative, so
%   it must never become an explicit diag. The result is
%   [2*(K+1) x 2*(K+1) x nDev], rows/cols ordered (re,im) within ascending
%   harmonic, which is the per-harmonic block ordering used here.

S = pagemtimes(GAMMA.*dRow,invGAMMA);

end


function rows = idxTermRow_ReIm(idxVar,numTerms,K)
%IDXTERMROW_REIM  Row indices of one variable inside a dense harmonic block.
%   The harmonic stride is 2*numTerms and (re,im) sit innermost, i.e.
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


function Xflat = hbDftPages(GAMMA,K,xTd)
%HBDFTPAGES  Harmonic coefficients of a page-based time-domain signal.
%   xTd is [nVar x N x nDev]; the result is [2*nVar*(K+1) x nDev], rows
%   harmonic-outer / variable-inner with (re,im) innermost. hbDft is not
%   page-aware, so the pages are folded into its variable dimension and the
%   whole fleet goes through one 2-D call. (The GENROU kernel's copy of this
%   helper is spelled dftPagesFlat; same code, same contract.)

[nVar,Nt,n] = size(xTd);
% Permuting to (var,page,time) keeps an inverter's variables contiguous in
% the fold, so splitting dim 1 afterwards recovers the page dimension
X2 = reshape(permute(xTd,[1 3 2]),nVar*n,Nt);
X_ReIm = hbDft(GAMMA,K,X2);             % [2*nVar*n x (K+1)]
X3 = reshape(X_ReIm,2*nVar,n,K+1);
% Put harmonics ahead of inverters, then flatten to harmonic-outer rows
Xflat = reshape(permute(X3,[1 3 2]),2*nVar*(K+1),n);

end
