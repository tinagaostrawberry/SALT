function [JvalsNlin,Fnlin,JvalsVabc,FVabc,dFdw0] = ...
    calcHbLoadCompMatrices_FD( ...
    ... % Time-domain previous guess (pages: [nVar x N x numLoadComp])
    vsTd,iphABCTd,wslipTd,deltaThetaSysTd,isTd,irTd,wrTd,PQTd,VmTd,izipTd, ...
    ... % Terminal bus voltage, time domain (pages: [3 x N x numLoadComp])
    vabcTd, ...
    ... % Frequency-domain previous guess, used for J*x
    w0, ...
    ... % Parameters (pages: [... x numLoadComp])
    L_PU,L_Te,T0,m, PQz,PQi,PQp, ...
    ... % Base changes
    busToLoadCompV,loadCompToBusI, ...
    ... % Transform
    K,N,GAMMA,invGAMMA,tSamp)
%CALCHBLOADCOMPMATRICES_FD  Nonlinear harmonic-balance composite-load kernel.
%
%   Dense-across-harmonics part of the composite load (induction motor +
%   ZIP): each derivative is formed in the time domain and squeezed back
%   into harmonic coefficients by the DFT sandwich
%   (GAMMA .* dRow) * invGAMMA. Every input carries the loads along the
%   third dimension, so the whole fleet is evaluated in one pass.
%   JvalsNlin is [nValsNlin x numLoadComp] with the dense blocks flattened
%   column-major and concatenated in the fixed order documented below;
%   Fnlin is the nonlinear residual [34*(K+1) x numLoadComp]. JvalsVabc is
%   the Park-V coupling into the terminal bus voltage, which is nonlinear
%   here and so cannot be handed back by the linear kernel.
%
%   Generalised from stampOfCompositeLoadMatrices_SS.m. SALT short-circuits
%   both Park transforms by setting theta_sys = 0 (vdq_s = vphA), which only
%   works at one harmonic; HB uses the genuine transform with
%   theta_sys(t) = w0*t + deltaThetaSys(t), which makes them DFT-sandwich
%   nonlinearities exactly like the GENROU's inverse Park.
%
%   PARK CONVENTION, AND WHY IT IS THE -sin ONE. With
%   theta_phAbc = theta + [0, -2*pi/3, +2*pi/3]:
%
%       T2    =  2/3*[ cos(theta_abc); -sin(theta_abc)]  (2x3, abc -> dq)
%       JT2   =  2/3*[-sin(theta_abc); -cos(theta_abc)]  = d(T2)/d(theta)
%       invT2 =      [ cos(theta_abc), -sin(theta_abc)]  (3x2, dq -> abc)
%       JinvT2=      [-sin(theta_abc), -cos(theta_abc)]  = d(invT2)/d(theta)
%
%   invT2 = (3/2)*T2.' is the exact inverse (T2*invT2 = I2), and each COLUMN
%   of invT2 sums to zero over the three phases, which is what makes
%   i_a+i_b+i_c = 0 hold by construction once the terminal current is
%   written abc FROM dq - SALT's +/-120 degree balance equations are then
%   unnecessary and are gone.
%
%   THE SIGN OF sin IS NOT FREE: it is fixed by the speed-voltage matrix the
%   rest of the model already carries. With this -sin convention and
%   theta_sys(t) = w0*t + deltaThetaSys(t), a balanced fundamental quantity
%   whose phase-A phasor is V maps to the dq pair
%
%       (vd,vq) = (Re{V*exp(-j*deltaThetaSys)}, Im{V*exp(-j*deltaThetaSys)})
%
%   i.e. the dq pair IS the phasor, rotated. The dq image of d/dt is then
%   +j*w, which is exactly SALT's wMat = [0 -w; w 0] and the ZIP's
%   J = [0 1; -1 0]; so every algebraic form in
%   stampOfCompositeLoadMatrices_SS carries over verbatim, which is the port
%   rule, and so do the initial values (SALT sets vdq_s = vphA and
%   thetaSys = 0, which this convention reproduces exactly).
%
%   Taking the OTHER convention, +sin as in getParkMatrices_2x3_PSCAD, makes
%   the dq pair the CONJUGATE of the phasor, so the dq image of d/dt becomes
%   -j*w and SALT's wMat, its ZIP J-matrix and its Q definition would all
%   have to flip sign - and the init values would need conjugating. That is
%   a self-consistent model too, and it is in fact the one SALT's GFM
%   inverter uses (its 2x2 invT = [cos sin; sin -cos] is precisely
%   "conjugate, then rotate", and its passive equations carry -stampReIm(jw)
%   where this device carries +wMat). But the composite load is written the
%   other way round, so mixing the two would silently reverse every
%   speed-voltage term in the machine.
%
%   See also STAMPOFHBLOADCOMPMATRICES_FD, GETINVPARKMATRICES,
%   CALCHBGENROUMATRICES_FD, CALCHBIBRGFMMATRICES_FD.

% VARIABLE LAYOUT (fixed): one composite-load block per harmonic is 34
% rows, every quantity a (real,imag) pair, so a global row is
% k*34 + localIndex, i.e. harmonic-outer / variable-inner, and a dq pair is
% stored (d_re, d_im, q_re, q_im). idxTermRow_ReIm below builds exactly that
% ordering inside each dense block.
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
% The stator/rotor equation is one 4x4 block over isr = [is; ir]
idxIsr = [idxHbLoadComp_Is, idxHbLoadComp_Ir];

% K is # harmonics, NOT including zero, so the total is +1
K_tot = K + 1;
% Loads run along dim 3
n = size(wrTd,3);
% Time runs along dim 2; the scalar states carry one row and tSamp is a
% column of samples (the page-based equivalent of the row-vector theta the
% 2-D Park helper asserts).
assert(size(wrTd,1)==1 && size(wslipTd,1)==1 && size(deltaThetaSysTd,1)==1)
assert(size(VmTd,1)==1 && size(tSamp,2)==1)
assert(size(vsTd,1)==2 && size(isTd,1)==2 && size(irTd,1)==2)
assert(size(PQTd,1)==2 && size(izipTd,1)==2)
assert(size(iphABCTd,1)==3 && size(vabcTd,1)==3)
assert(size(vsTd,2)==N && size(vabcTd,2)==N && size(tSamp,1)==N)
assert(all(cellfun(@(x)size(x,3)==n, ...
    {vsTd,iphABCTd,wslipTd,deltaThetaSysTd,isTd,irTd,PQTd,VmTd,izipTd, ...
     vabcTd,L_PU,L_Te,T0,m,PQz,PQi,PQp,busToLoadCompV,loadCompToBusI})))
% Each DFT sandwich must come out [2*K_tot x 2*K_tot] to land on a harmonic
% block, which holds only when GAMMA keeps exactly DC..K (N = 2*K+1).
assert(size(GAMMA,1)==2*K_tot && size(GAMMA,2)==N)
assert(size(invGAMMA,1)==N && size(invGAMMA,2)==2*K_tot)

% Number of (re,im) pairs addressed by each group, i.e. the "numTerms" that
% idxTermRow_ReIm needs to work out its harmonic stride
nTermsDq   = numel(idxHbLoadComp_Vs)/2;             % 2: a dq pair
nTermsAbc  = numel(idxHbLoadComp_IphABC)/2;         % 3: three phases
nTermsIsr  = numel(idxIsr)/2;                       % 4: (d,q) stator+rotor
nTermsScal = 1;                                     % a scalar quantity

% Block layout, in the order they are concatenated into JvalsNlin. Each is a
% dense (numel(rows)*K_tot) x (numel(cols)*K_tot) array.
%  # | rows  | cols  | term
%  1 | 1:4   | 13:14 | d(T2(th)*v_abc)/d(deltaThetaSys)      ->  8*K_tot^2
%  2 | 5:10  | 15:18 | d(invT2(th)*(is+izip))/d(is)          -> 24*K_tot^2
%  3 | 5:10  | 31:34 | d(invT2(th)*(is+izip))/d(izip)        -> 24*K_tot^2
%  4 | 5:10  | 13:14 | d(invT2(th)*(is+izip))/d(deltaTheta)  -> 12*K_tot^2
%  5 | 15:22 | 15:22 | d(wMatSlip*L_PU*isr)/d(isr)           -> 64*K_tot^2
%  6 | 15:22 | 11:12 | d(wMatSlip*L_PU*isr)/d(wslip)         -> 16*K_tot^2
%  7 | 23:24 | 15:22 | d(isr'*L_Te*isr)/d(isr)               -> 16*K_tot^2
%  8 | 23:24 | 23:24 | d(-T0*(wr/w0)^m)/d(wr)                ->  4*K_tot^2
%  9 | 25:28 | 29:30 | d(PQz*Vm^2 + PQi*Vm)/d(Vm)            ->  8*K_tot^2
% 10 | 29:30 | 1:4   | d(|vs|)/d(vs)                         ->  8*K_tot^2
% 11 | 31:34 | 25:28 | d(iZipEval/Vm^2)/d(PQ)                -> 16*K_tot^2
% 12 | 31:34 | 1:4   | d(iZipEval/Vm^2)/d(vs)                -> 16*K_tot^2
% 13 | 31:34 | 29:30 | d(iZipEval/Vm^2)/d(Vm)                ->  8*K_tot^2
% Sizes are taken from the index vectors so the offsets cannot drift from
% the documented layout.
nBlkVals = K_tot^2 * [ ...
    numel(idxHbLoadComp_Vs)*numel(idxHbLoadComp_DeltaThetaSys); ...     %  1
    numel(idxHbLoadComp_IphABC)*numel(idxHbLoadComp_Is); ...            %  2
    numel(idxHbLoadComp_IphABC)*numel(idxHbLoadComp_Izip); ...          %  3
    numel(idxHbLoadComp_IphABC)*numel(idxHbLoadComp_DeltaThetaSys); ... %  4
    numel(idxIsr)^2; ...                                               %  5
    numel(idxIsr)*numel(idxHbLoadComp_WSlip); ...                       %  6
    numel(idxHbLoadComp_Wr)*numel(idxIsr); ...                          %  7
    numel(idxHbLoadComp_Wr)^2; ...                                      %  8
    numel(idxHbLoadComp_PQ)*numel(idxHbLoadComp_Vm); ...                %  9
    numel(idxHbLoadComp_Vm)*numel(idxHbLoadComp_Vs); ...                % 10
    numel(idxHbLoadComp_Izip)*numel(idxHbLoadComp_PQ); ...              % 11
    numel(idxHbLoadComp_Izip)*numel(idxHbLoadComp_Vs); ...              % 12
    numel(idxHbLoadComp_Izip)*numel(idxHbLoadComp_Vm)];                 % 13
assert(isequal(nBlkVals.', ...
    K_tot^2*[8 24 24 12 64 16 16 4 8 8 16 16 8]))
% Running start offset of each block inside a JvalsNlin column
offBlk = [0; cumsum(nBlkVals)];
nValsNlin = offBlk(end);
assert(nValsNlin==224*K_tot^2)
JvalsNlin = zeros(nValsNlin,n);

% ----------------------- Park transform, in time -------------------------
% theta_sys(t) = w0*t + deltaThetaSys(t). tSamp is transposed so time lines
% up with dim 2 of the page arrays.
thetaTd = w0*(tSamp.') + deltaThetaSysTd;                   % [1 x N x n]
% The three phase angles are kept as explicit [3 x N x n] rows rather than
% concatenated along time (as the 2-D Park helpers do), because a row vector
% of angles cannot carry a page dimension.
thetaTd_phAbc = thetaTd + [0;-2*pi/3;2*pi/3];               % [3 x N x n]
cosTd = cos(thetaTd_phAbc);
sinTd = sin(thetaTd_phAbc);
% Forward Park rows (2x3) and their theta-derivative:
%   T2 = 2/3*[cos(theta_abc); -sin(theta_abc)],  JT2 = d(T2)/d(theta)
% Stored as a cell of ROWS so T2row{i}(p,:,:) is the (i,p) entry over time.
T2row  = {(2/3)*cosTd, -(2/3)*sinTd};
JT2row = {-(2/3)*sinTd, -(2/3)*cosTd};
% Inverse Park columns (3x2) = (3/2)*T2.' = [cos, -sin], and its
% theta-derivative [-sin, -cos]. Stored as a cell of COLUMNS so
% invT2col{j}(p,:,:) is the (p,j) entry over time.
invT2col  = {cosTd, -sinTd};
JinvT2col = {-sinTd, -cosTd};

% ------------------------------ Park V -----------------------------------
% -vs + T2(thetaSys)*v_abc*busToLoadCompV = 0
% (Val) row i is the i-th T2 row dotted with v_abc at each time sample
ParkV = cat(1, ...
    sum(T2row{1}.*vabcTd,1), ...
    sum(T2row{2}.*vabcTd,1)).*busToLoadCompV;               % [2 x N x n]
% (Deriv) w.r.t. deltaThetaSys: d(theta_sys)/d(deltaThetaSys) = 1, so the
% chain rule is just JT2 in place of T2.
dParkV_dDeltaTheta = zeros(2*nTermsDq*K_tot,2*nTermsScal*K_tot,n);
for idxVar = 1:nTermsDq
    dParkV_dDeltaTheta(idxTermRow_ReIm(idxVar,nTermsDq,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA, ...
        sum(JT2row{idxVar}.*vabcTd,1).*busToLoadCompV);
end
% (Deriv) w.r.t. the three bus-voltage nodes. This is the coupling OUT of
% the device block, so it is returned separately. Its columns are ordered
% harmonic-outer / node-inner, the same ordering idxTermRow_ReIm produces
% for three terms, which is how the caller lays the bus columns out.
nValsVabc = numel(idxHbLoadComp_Vs)*2*nTermsAbc*K_tot^2;
assert(nValsVabc==24*K_tot^2)
dParkV_dVabc = zeros(2*nTermsDq*K_tot,2*nTermsAbc*K_tot,n);
for idxVar = 1:nTermsDq
    for jdxVar = 1:nTermsAbc
        dParkV_dVabc(idxTermRow_ReIm(idxVar,nTermsDq,K), ...
            idxTermRow_ReIm(jdxVar,nTermsAbc,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            T2row{idxVar}(jdxVar,:,:).*busToLoadCompV);
    end
end
JvalsVabc = reshape(dParkV_dVabc,nValsVabc,n);

% ------------------------------ Park I -----------------------------------
% -i_abc + invT2(thetaSys)*(is + izip)*loadCompToBusI = 0
% The motor and the ZIP share one terminal, so the inverse Park acts on
% their SUM; the base change runs load -> bus because iphABC is in BUS
% per-unit while the dq currents are in LOAD per-unit.
isizTd = isTd + izipTd;                                     % [2 x N x n]
% (Val) phase p is row p of invT2 dotted with (is+izip) at each time sample
ParkI = (invT2col{1}.*isizTd(1,:,:) + ...
         invT2col{2}.*isizTd(2,:,:)).*loadCompToBusI;       % [3 x N x n]
% (Deriv) w.r.t. is. The matrix is the same for izip, so the block is
% computed once and reused - exactly as the GENROU reuses its inverse-Park
% voltage derivative for the current.
dParkI_dIs = zeros(2*nTermsAbc*K_tot,2*nTermsDq*K_tot,n);
for idxVar = 1:nTermsAbc
    for jdxVar = 1:nTermsDq
        dParkI_dIs(idxTermRow_ReIm(idxVar,nTermsAbc,K), ...
            idxTermRow_ReIm(jdxVar,nTermsDq,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            invT2col{jdxVar}(idxVar,:,:).*loadCompToBusI);
    end
end
dParkI_dIzip = dParkI_dIs;
% (Deriv) w.r.t. deltaThetaSys: d(invT2*x)/d(theta) contracted over dq at
% each time sample
dParkI_dDeltaTheta = zeros(2*nTermsAbc*K_tot,2*nTermsScal*K_tot,n);
for idxVar = 1:nTermsAbc
    dParkI_dDeltaTheta(idxTermRow_ReIm(idxVar,nTermsAbc,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA, ...
        (JinvT2col{1}(idxVar,:,:).*isizTd(1,:,:) + ...
         JinvT2col{2}(idxVar,:,:).*isizTd(2,:,:)).*loadCompToBusI);
end

% ------------------- Stator/rotor speed voltage --------------------------
% The stator half of wMat uses wSys == w0, a single DC scalar, so it is
% linear and harmonic-diagonal and lives in the linear kernel. Only the
% ROTOR half is nonlinear, because wslip is harmonic-expanded:
%   wMatSlip*L_PU*isr,  wMatSlip(3,4) = -wslip, wMatSlip(4,3) = +wslip
% i.e. it touches only the two rotor equations, rows 19:22.
isrTd   = [isTd; irTd];                                     % [4 x N x n]
phisrTd = pagemtimes(L_PU,isrTd);                           % flux linkage
% Term rows 3 and 4 of the 4-row block are the rotor pair; each takes MINUS
% wslip times the OTHER rotor flux (row 3) or PLUS it (row 4), which is
% SALT's wMat(3,4) = -wSlip, wMat(4,3) = +wSlip.
rotTerm = [3 4];
rotPrt  = [4 3];                % the flux row each rotor equation picks up
rotSgn  = [-1 1];
% (Val)
SlipSpeedVolt = zeros(4,N,n);
for q = 1:2
    SlipSpeedVolt(rotTerm(q),:,:) = ...
        rotSgn(q)*wslipTd.*phisrTd(rotPrt(q),:,:);
end
% (Deriv) w.r.t. isr: the flux row is linear in isr, so the derivative is
% wslip times the corresponding L_PU entry, held constant over time.
dSlipSr_dIsr = zeros(2*nTermsIsr*K_tot,2*nTermsIsr*K_tot,n);
for q = 1:2
    for jdxVar = 1:nTermsIsr
        dSlipSr_dIsr(idxTermRow_ReIm(rotTerm(q),nTermsIsr,K), ...
            idxTermRow_ReIm(jdxVar,nTermsIsr,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            rotSgn(q)*wslipTd.*L_PU(rotPrt(q),jdxVar,:));
    end
end
% (Deriv) w.r.t. wslip: SALT's [zeros(2,1,n); -phisr(4); phisr(3)]
dSlipSr_dWslip = zeros(2*nTermsIsr*K_tot,2*nTermsScal*K_tot,n);
for q = 1:2
    dSlipSr_dWslip(idxTermRow_ReIm(rotTerm(q),nTermsIsr,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA, ...
        rotSgn(q)*phisrTd(rotPrt(q),:,:));
end

% ------------------------- Swing equation --------------------------------
% isr'*L_Te*isr - T0*(wr/wSys)^m, with wSys the scalar w0. Both terms are
% nonlinear; the restored inertia term -(2H/w0nom)*j*k*w0*wr is linear and
% lives in the linear kernel.
% (Deriv) of the quadratic form: d/d(isr) (isr'*L_Te*isr) = (L_Te+L_Te')*isr
% at each time sample, SALT's LTeSym.
dte_disr = pagemtimes(L_Te + pagetranspose(L_Te),isrTd);    % [4 x N x n]
dSwing_dIsr = zeros(2*nTermsScal*K_tot,2*nTermsIsr*K_tot,n);
for jdxVar = 1:nTermsIsr
    dSwing_dIsr(:,idxTermRow_ReIm(jdxVar,nTermsIsr,K),:) = ...
        dftSandwich(GAMMA,invGAMMA,dte_disr(jdxVar,:,:));
end
% (Deriv) of the load torque, SALT's -T0*(m*wr^(m-1))/(wSys^m) with
% wSys = w0
dSwing_dWr = dftSandwich(GAMMA,invGAMMA, ...
    -T0.*(m.*wrTd.^(m-1))./(w0.^m));
% (Val) the quadratic form is per time sample, so contract dim 1 rather than
% calling pagemtimes with a time-indexed middle dimension
loadTorqueTd = T0.*(wrTd./w0).^m;
Swing = sum(isrTd.*pagemtimes(L_Te,isrTd),1) - loadTorqueTd;
% (Deriv) of the load torque w.r.t. w0. THIS BELONGS HERE, NOT IN THE LINEAR
% KERNEL, and the reason is the N scaling. The swing residual carries the
% load torque with a MINUS, and pointwise in time
%
%   d/dw0[-T0*(wr(t)/w0)^m] = +T0*m*(wr(t)/w0)^m/w0 = +(m/w0)*loadTorque(t)
%
% - PLUS, which is SALT's dFswing_dwSys = -T0*wr^m*(-m)*wSys^(-m-1)
% rearranged. So the exact column is +(m/w0) times the DFT of the load
% torque itself, a quantity only this kernel has, because it is the one that
% works in the time domain. Writing the same algebra on the harmonic
% COEFFICIENT instead, as the linear kernel first did, is wrong by N^(m-1):
% a DC coefficient is N times the time value, the m-th power multiplies that
% error m times over, and the single 1/w0 divides it out only once. With
% m = 2 and N = 3 that is a factor of 3 on this one Jacobian entry, which
% was enough to turn Newton's quadratic convergence into a linear crawl -
% and, because this column also feeds the right-hand side, to move the
% converged answer.
% m is a page, the DFT result a [2*K_tot x n] column set, so m is flattened
% to a row first and broadcasts down the harmonics.
dSwing_dw0 = (reshape(m,1,n)./w0).*dftPagesFlat(GAMMA,K,loadTorqueTd);

% ---------------------------- ZIP power ----------------------------------
% -PQ + PQz*Vm^2 + PQi*Vm + PQp = 0
% (Val) PQp is a pure CONSTANT and is carried here rather than in the linear
% kernel's F: run through the DFT it lands in the DC bin already multiplied
% by N, which is the scaling the unnormalized coefficient state vector
% needs, and the caller's J*x - Fnlin then puts it on the right-hand side
% with the right sign.
PQzipVal = PQz.*VmTd.^2 + PQi.*VmTd + PQp;                  % [2 x N x n]
% (Deriv) SALT's 2*PQz*Vm + PQi
dPQzip_dVm_td = 2*PQz.*VmTd + PQi;
dZip_dVm = zeros(2*nTermsDq*K_tot,2*nTermsScal*K_tot,n);
for idxVar = 1:nTermsDq
    dZip_dVm(idxTermRow_ReIm(idxVar,nTermsDq,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA,dPQzip_dVm_td(idxVar,:,:));
end

% ------------------------ ZIP voltage magnitude --------------------------
% -(1 + j*k*w0*Tau)*Vm + |vs| = 0; only |vs| is nonlinear.
% SALT: Vm = (vd_s^2+vq_s^2)^(1/2)
vsMagTd = sqrt(vsTd(1,:,:).^2 + vsTd(2,:,:).^2);            % [1 x N x n]
% (Deriv) SALT's pagetranspose(vs)./vsMag, i.e. [vd vq]/|vs|
dVsMag_dVs_td = vsTd./vsMagTd;                              % [2 x N x n]
dVsMag_dVs = zeros(2*nTermsScal*K_tot,2*nTermsDq*K_tot,n);
for jdxVar = 1:nTermsDq
    dVsMag_dVs(:,idxTermRow_ReIm(jdxVar,nTermsDq,K),:) = ...
        dftSandwich(GAMMA,invGAMMA,dVsMag_dVs_td(jdxVar,:,:));
end

% ---------------------------- ZIP current --------------------------------
% SALT: idZip = (P*vd + Q*vq)/Vm^2
%       iqZip = (P*vq - Q*vd)/Vm^2
% i.e. izip = [vs, J*vs]*PQ / Vm^2 with J = [0 1; -1 0]
iZipEval = cat(1, ...
    vsTd(1,:,:).*PQTd(1,:,:) + vsTd(2,:,:).*PQTd(2,:,:), ...
    vsTd(2,:,:).*PQTd(1,:,:) - vsTd(1,:,:).*PQTd(2,:,:));   % [2 x N x n]
Vm2_inv  = 1./VmTd.^2;
IzipVal  = iZipEval.*Vm2_inv;
% (Deriv) SALT's blk, split into its three column groups. Rows are the
% (id,iq) equations, so entry (i,j) is read off
%   blk = [ vd  vq |  P  Q ]
%         [ vq -vd | -Q  P ]  all times 1/Vm^2
% and the Vm column is -2*iZipEval/Vm^3.
dIzip_dPQ_td = {{ vsTd(1,:,:), vsTd(2,:,:)}, ...            % row id
                { vsTd(2,:,:),-vsTd(1,:,:)}};               % row iq
dIzip_dVs_td = {{ PQTd(1,:,:), PQTd(2,:,:)}, ...            % row id
                {-PQTd(2,:,:), PQTd(1,:,:)}};               % row iq
dIzip_dVm_td = -2*iZipEval./VmTd.*Vm2_inv;
dIzip_dPQ = zeros(2*nTermsDq*K_tot,2*nTermsDq*K_tot,n);
dIzip_dVs = zeros(2*nTermsDq*K_tot,2*nTermsDq*K_tot,n);
dIzip_dVm = zeros(2*nTermsDq*K_tot,2*nTermsScal*K_tot,n);
for idxVar = 1:nTermsDq
    for jdxVar = 1:nTermsDq
        dIzip_dPQ(idxTermRow_ReIm(idxVar,nTermsDq,K), ...
            idxTermRow_ReIm(jdxVar,nTermsDq,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            dIzip_dPQ_td{idxVar}{jdxVar}.*Vm2_inv);
        dIzip_dVs(idxTermRow_ReIm(idxVar,nTermsDq,K), ...
            idxTermRow_ReIm(jdxVar,nTermsDq,K),:) = ...
            dftSandwich(GAMMA,invGAMMA, ...
            dIzip_dVs_td{idxVar}{jdxVar}.*Vm2_inv);
    end
    dIzip_dVm(idxTermRow_ReIm(idxVar,nTermsDq,K),:,:) = ...
        dftSandwich(GAMMA,invGAMMA,dIzip_dVm_td(idxVar,:,:));
end

% ------------------------ Jacobian values --------------------------------
JvalsNlin(offBlk( 1)+(1:nBlkVals( 1)),:) = ...
    reshape(dParkV_dDeltaTheta,nBlkVals( 1),n);
JvalsNlin(offBlk( 2)+(1:nBlkVals( 2)),:) = ...
    reshape(dParkI_dIs,nBlkVals( 2),n);
JvalsNlin(offBlk( 3)+(1:nBlkVals( 3)),:) = ...
    reshape(dParkI_dIzip,nBlkVals( 3),n);
JvalsNlin(offBlk( 4)+(1:nBlkVals( 4)),:) = ...
    reshape(dParkI_dDeltaTheta,nBlkVals( 4),n);
JvalsNlin(offBlk( 5)+(1:nBlkVals( 5)),:) = ...
    reshape(dSlipSr_dIsr,nBlkVals( 5),n);
JvalsNlin(offBlk( 6)+(1:nBlkVals( 6)),:) = ...
    reshape(dSlipSr_dWslip,nBlkVals( 6),n);
JvalsNlin(offBlk( 7)+(1:nBlkVals( 7)),:) = ...
    reshape(dSwing_dIsr,nBlkVals( 7),n);
JvalsNlin(offBlk( 8)+(1:nBlkVals( 8)),:) = ...
    reshape(dSwing_dWr,nBlkVals( 8),n);
JvalsNlin(offBlk( 9)+(1:nBlkVals( 9)),:) = ...
    reshape(dZip_dVm,nBlkVals( 9),n);
JvalsNlin(offBlk(10)+(1:nBlkVals(10)),:) = ...
    reshape(dVsMag_dVs,nBlkVals(10),n);
JvalsNlin(offBlk(11)+(1:nBlkVals(11)),:) = ...
    reshape(dIzip_dPQ,nBlkVals(11),n);
JvalsNlin(offBlk(12)+(1:nBlkVals(12)),:) = ...
    reshape(dIzip_dVs,nBlkVals(12),n);
JvalsNlin(offBlk(13)+(1:nBlkVals(13)),:) = ...
    reshape(dIzip_dVm,nBlkVals(13),n);

% ----------------------- Nonlinear residual ------------------------------
% Fnlin holds the harmonic coefficients of the nonlinear terms themselves,
% g(x_prev), at the rows the terms belong to and zero elsewhere. With the
% caller's J*x - F convention the caller forms F = JvalsNlin*x_prev - Fnlin,
% so the minus on g(x_prev) lives in the caller and Fnlin is returned with
% the plain sign of the term - the same contract as the GENROU kernel.
Fnlin = zeros(nEqns*K_tot,n);
% Rows are gathered harmonic-outer / variable-inner, matching the layout
% that the DFT sandwich and hbDft both produce.
Fnlin(hbRows(idxHbLoadComp_Vs,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,ParkV);
Fnlin(hbRows(idxHbLoadComp_IphABC,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,ParkI);
Fnlin(hbRows(idxIsr,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,SlipSpeedVolt);
Fnlin(hbRows(idxHbLoadComp_Wr,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,Swing);
Fnlin(hbRows(idxHbLoadComp_PQ,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,PQzipVal);
Fnlin(hbRows(idxHbLoadComp_Vm,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,vsMagTd);
Fnlin(hbRows(idxHbLoadComp_Izip,K,nEqns),:) = ...
    dftPagesFlat(GAMMA,K,IzipVal);
% Rows 11:12 (slip) and 13:14 (frame angle) stay zero: both equations are
% entirely linear, so they have no nonlinear residual.

% The Park V residual is already in Fnlin rows 1:4 - JvalsVabc is only the
% Jacobian column group that leaves the device block, and the residual is
% not split along with it. So there is nothing for the caller to add here.
FVabc = [];

% ----------------- Coupling into the system-frequency column -------------
% The part of the w0 column that can only be formed in the time domain:
% the load torque's. Everything else in that column is algebraic in the
% harmonic coefficients and stays in the linear kernel. Same shape and
% layout as Fnlin, so the caller just adds the two columns together; and
% like the linear kernel's dFdw0, this half IS state-dependent, so its
% dFdw0*w0 goes into the right-hand side.
dFdw0 = zeros(nEqns*K_tot,n);
dFdw0(hbRows(idxHbLoadComp_Wr,K,nEqns),:) = dSwing_dw0;

% iphABCTd is unused: the terminal current enters the Park I equation only
% through the linear -i_abc term, so it has no nonlinear derivative. It is
% taken as an argument so both kernels share one state list, and so its
% shape and page count are checked above.

end


function S = dftSandwich(GAMMA,invGAMMA,dRow)
%DFTSANDWICH  (GAMMA .* dRow) * invGAMMA for every load at once.
%   dRow is [1 x N x numLoadComp]: one time-domain derivative row per load.
%   The broadcast ".*" IS the diagonal multiplication by the derivative, so
%   it must never become an explicit diag. The result is
%   [2*(K+1) x 2*(K+1) x numLoadComp], rows/cols ordered (re,im) within
%   ascending harmonic, which is the per-harmonic block ordering used here.

S = pagemtimes(GAMMA.*dRow,invGAMMA);

end


function rows = idxTermRow_ReIm(idxVar,numTerms,K)
%IDXTERMROW_REIM  Row indices of one variable inside a dense harmonic block.
%   The harmonic stride is 2*numTerms and (re,im) sit innermost, i.e.
%   harmonic-outer / variable-inner - the same ordering the device rows use.

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
%   xTd is [nVar x N x numLoadComp]; the result is
%   [2*nVar*(K+1) x numLoadComp], rows harmonic-outer / variable-inner with
%   (re,im) innermost. hbDft is NOT page-aware, so the pages are folded into
%   its variable dimension and unfolded afterwards.

[nVar,Nt,n] = size(xTd);
% GAMMA acts only along time, so the pages fold into hbDft's variable
% dimension and the whole fleet goes through one 2-D call. Permuting to
% (var,page,time) keeps a load's variables contiguous in the fold.
X2 = reshape(permute(xTd,[1 3 2]),nVar*n,Nt);
X_ReIm = hbDft(GAMMA,K,X2);             % [2*nVar*n x (K+1)]
% Undo the fold: within a column the rows run (re,im) inside variable, then
% load, so splitting dim 1 recovers the page dimension.
X3 = reshape(X_ReIm,2*nVar,n,K+1);
% Put harmonics ahead of loads, then flatten to harmonic-outer rows
Xflat = reshape(permute(X3,[1 3 2]),2*nVar*(K+1),n);

end
