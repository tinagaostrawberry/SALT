classdef HB < handle
%HB  Harmonic-balance steady-state solver for three-phase EMT models.
%
%   hb = HB(baseMVA,bus,slack_busNum,slack_imScalConst,numHarmonics,maxDeltaT)
%
%   HB finds the PERIODIC steady state of the full three-phase EMT
%   differential-algebraic system,
%
%       x(t) = sum_{k=0..K} X_k exp(j*k*w0*t) + conj(...)
%
%   with the fundamental w0 itself an unknown. That makes it an AUTONOMOUS
%   harmonic-balance problem - a free-running oscillator rather than a
%   driven one - which is what the extra slack-angle equation is for
%   (see STAMPPHASE_Y_I).
%
%   Relationship to SALT. SALT is the K = 1, balanced case: at a single
%   harmonic the dq frame is stationary, so every dq quantity is a real DC
%   scalar and the machine reduces to (id,iq,ifd) with the damper currents
%   identically zero. HB cannot make that reduction - harmonic flux induces
%   damper currents, and the dampers are the dominant harmonic path through
%   a machine - so it carries all seven windings with a (real,imag) pair at
%   every harmonic. With no harmonic excitation and K = 1 the harmonic-1 dq
%   components go to zero, the dampers with them, and HB's DC components
%   reduce to SALT's answer.
%
%   Method. Three layers, which is the whole architecture:
%
%     linear network    Exact per-harmonic phasor admittances, evaluated at
%                       w_k = k*w0 in closed form. The network is LTI, so
%                       harmonic k sees only harmonic k and this part of
%                       the Jacobian is block diagonal in harmonic. The
%                       only off-block entries are the w0 column, which
%                       every harmonic pulls on with weight k.
%
%     phase gauge       A periodic solution of an autonomous system is
%                       invariant under time shift, so the Jacobian is
%                       singular by one dimension unless the phase is
%                       pinned. One equation fixes the slack bus's
%                       fundamental ANGLE and leaves its magnitude free.
%
%     nonlinear devices GAMMA*diag(dg/dx)*invGAMMA - evaluate the
%                       nonlinearity and its derivative POINTWISE IN TIME,
%                       then transform back. These blocks are dense in
%                       harmonic index, and that harmonic mixing is the
%                       point of the whole exercise.
%
%   Assembly follows SALT exactly: BUILDSTAMPMAPS precomputes the row and
%   column index of every Jacobian entry once, the stamping methods only
%   refill obj.stampVal, and ASSEMBLE_Y_I turns the three arrays into a
%   sparse matrix with duplicates summed.
%
%   ---------------------------------------------------------------------
%   HOW A SALT DEVICE MODEL BECOMES AN HB ONE
%   ---------------------------------------------------------------------
%   SALT solves one harmonic, so in its rotating frame every dq quantity
%   is constant and it DROPS every d/dt term. HB restores them. The
%   mechanical rule is:
%
%       take SALT's algebraic form verbatim, then add back each dropped
%       derivative as (j*k*w0) acting on the harmonic coefficient, keeping
%       whatever constant already stood in front of d/dt.
%
%   At k = 0 every restored term vanishes and HB reduces to SALT exactly,
%   which is what validation/validateHbVsSalt.m asserts.
%
%   THE TIME BASE IS NOT THE SAME IN EVERY BLOCK, and getting it wrong
%   scales a harmonic impedance by w0nom (~377) without any error. The
%   authority is the TD model, where the derivative is actually
%   discretised: a block written against dT is in physical seconds, one
%   written against dTPU = dT*w0nom is in normalised time. Per block:
%
%     GENROU stator        L_PU = L/w0nom, physical seconds
%                          -> j*k*w0*L_PU          (= j*k*L)
%     GENROU swing         2H/w0nom, physical seconds
%                          -> j*k*w0*2H/w0nom
%     GENROU angle         physical seconds        -> j*k*w0
%     EXDC1, TGOV1D        time constants in seconds
%                          -> j*k*w0*T
%     IM stator/rotor      L_PU = L/w0nom, physical seconds
%                          -> j*k*w0*L_PU
%     IM swing             2H/w0nom, physical seconds
%                          -> j*k*w0*2H/w0nom
%     IM ZIP magnitude     Tau in seconds; the SS kernel has NO filter
%                          because d/dt = 0 makes Vm = |vs|. HB restores
%                          it as (1 + j*k*w0*Tau)*Vm = |vs|, so the ZIP
%                          sees an attenuated voltage at every harmonic
%                          and behaves as an impedance there. With
%                          Tau = 0.1 s the 2nd-harmonic attenuation is
%                          about 1.3%.
%     IBR LCL + coupling   dTPU, i.e. NORMALISED time
%                          -> j*k*w0/w0nom
%     IBR inner loops      kc_p, kv_p against dT, physical seconds
%                          -> j*k*w0
%     IBR droop filters,   physical seconds        -> j*k*w0
%     theta, virtual
%     resistor filter
%
%   Nonlinear products go through the DFT sandwich rather than being
%   evaluated once at the operating point. Which quantities are
%   harmonic-expanded and which stay DC is a modelling commitment, not a
%   derivation: the system fundamental w0 is a SINGLE SCALAR, while every
%   rotor speed is expanded. So the induction motor's slip is
%
%       wSlip_k = w0*delta_k0 - wr_k
%
%   which is exactly how the GENROU already treats d(deltaTheta)/dt =
%   wr - w0, where w0 couples only at k = 0 ("no w0 at higher harmonics by
%   definition of the fundamental"). The stator speed-voltage term is then
%   bilinear in (w0, i) and needs a d/dw0 stamp, like a network branch.
%
%   A PRE-EXISTING CONVENTION TO BE AWARE OF: SALT's IBR divides the two
%   measurement bandwidths by 2*pi (w_VR = w_VR_rad/(2*pi), likewise
%   wmeas) in BOTH its steady-state and time-domain kernels, so it is
%   consistent rather than a transcription slip - but it means a data
%   value of 120 becomes a ~19 rad/s filter. HB reproduces it unchanged.
%   It matters more here than in SALT, because in HB those filters
%   actually roll something off.

    properties (Constant)
        %% System
        f0nom = 60                    % nominal system frequency, Hz
        w0nom = 2*pi*HB.f0nom
        % Spelled w0 as well, because SALT's composite-load and GFM-IBR
        % initialisation routines read obj.w0 and HB reuses them verbatim
        % rather than keeping a second, drifting copy of 500 lines of
        % operating-point arithmetic. Same value, nominal in both.
        w0    = 2*pi*HB.f0nom
        genP  = 2                     % two-pole machine; see TGOV1D kernel

        %% Equation counts, per device PER HARMONIC
        % Every quantity is carried as a (real,imag) pair, so each counts 2.
        numReImEqns     = 2
        numGenrouEqns   = 34          % 17 variables (see HBINITINDICES)
        numExcDc1Eqns   = 12          % 5 states + 1 vt calculation
        numTgov1dEqns   = 6           % 2 states + 1 pmech calculation
        % SALT carries 20 and 33 for these two; HB needs a (real,imag) pair
        % per harmonic for every dq quantity, so a dq pair costs 4 slots
        % where SALT spent 2. See HBINITINDICES.
        numLoadCompEqns = 34          % IM: 4 stator volt + 6 abc current +
                                      % 2 slip + 2 frame angle + 4 stator
                                      % current + 4 rotor current + 2 speed
                                      % ZIP: 4 power + 2 magnitude + 4 current
        numIbrGFMEqns   = 62          % 6 abc current + 8 dq pairs of the
                                      % passive network and controllers +
                                      % 2 speed + 2 Q droop + 4 measured PQ
                                      % + 2 frame angle + 12 virtual resistor
                                      % + 2 dq0 zero sequence (io0 and the
                                      % v0 = 0 row); see the DELTA, NOT WYE
                                      % note in CALCHBIBRGFMMATRICES_FD

        %% Sanity-check tolerance
        absTolCheck = 1e-9
    end

    properties
        %% General
        baseMVA
        bus
        numBus = 0
        caseName
        bus_BaseKv
        bus_IBase

        % Slack. The frequency unknown needs one extra equation, which pins
        % the slack angle, so store the imag/real ratio of its voltage.
        slack_busNum
        slack_imScalConst

        % Y*x_v = I, stamped as J(x[k])*x[k+1] = J(x[k])*x[k] - F(x[k])
        Y
        I
        x_v
        x_v0

        %% Harmonic-balance grid
        K                             % number of harmonics above DC
        f0                            % current fundamental, Hz
        N                             % DFT length
        tSamp
        fSamp
        fGamma
        GAMMA
        invGAMMA
        maxDeltaT
        tSim
        deltaT

        %% Network - transmission lines (three-phase PI, no mutual coupling)
        tx_busFrom
        tx_busTo
        tx_R
        tx_L
        tx_C
        tx_mnaIdxMap = zeros(0,2)
        tx_numMnaVars = 0

        %% Network - transformers (three-phase R-L, no mutual coupling)
        xfmr_busFrom
        xfmr_busTo
        xfmr_R
        xfmr_L
        xfmr_mnaIdxMap = zeros(0,2)
        xfmr_numMnaVars = 0

        %% Network - constant impedances to ground
        load_bus
        load_R
        load_X
        load_Z
        load_mnaIdxMap = zeros(0,2)
        load_numMnaVars = 0

        %% Harmonic current injection
        % Without an excitation source the harmonics decouple from the
        % fundamental and the only solution is the fundamental-only one, so
        % something has to drive them.
        includeHigherHarmonics = false
        higherHarmonics_bus
        % NOTE: the legacy called this "percentFundamental" but used it as a
        % bare multiplier, so the value 2 means 200% of the power-flow
        % fundamental current, not 2%. Renamed to say what it does; the
        % arithmetic is unchanged.
        higherHarmonics_scaleFundamental
        higherHarmonics_harmonicNum
        higherHarmonics_I

        %% GENROU
        numGenrou = 0
        genrou_bus
        genrou_EsBase
        genrou_IsBase
        genrou_IsBaseRMS
        genrou_ZsBase
        genrou_LsBase
        genrou_BaseMVA
        genrou_MW0
        genrou_Mvar0
        genrou_RMatrix                % [7 x 7 x numGenrou]
        genrou_LMatrix
        genrou_LspMatrix
        genrou_H
        genrou_KD
        genrou_TorqueM0
        genrou_efd0
        genrou_Lad
        genrou_Rfd
        genrou_Xd
        genrou_Xq
        genrou_Ra

        %% EXDC1
        % All controller parameters are assumed to be on the MACHINE base,
        % not the system base. See the note in SALT; the IEEE PES governor
        % report says the true base is the turbine rating, which is not in
        % the data, so the machine rating stands in for it.
        numExcDc1 = 0
        excDc1_bus
        excDc1_Tr
        excDc1_Ta
        excDc1_Tc
        excDc1_Tb
        excDc1_Te
        excDc1_Tf1
        excDc1_Kf1
        excDc1_Ka
        excDc1_Ke
        excDc1_vref0
        excDc1_genToExc               % Lad/Rfd, the field-circuit gain
        excDc1_genrouIdx

        %% TGOV1D
        numTgov1d = 0
        tgov1d_bus
        tgov1d_T1
        tgov1d_T2
        tgov1d_T3
        tgov1d_R
        tgov1d_Dt
        tgov1d_Pref0
        tgov1d_genrouIdx

        %% COMPOSITE LOAD (induction motor + ZIP)
        % Property NAMES mirror SALT's exactly. That is load-bearing, not
        % cosmetic: HB reuses SALT's compLoadValInitCalculations_assumeW0
        % and ibrGFMValInitCalculations_assumeW0 verbatim, and those read
        % obj.<name> directly.
        numLoadComp = 0
        loadComp_bus
        loadComp_VBase
        loadComp_IBase
        loadComp_IsBaseRMS
        loadComp_MVA_base
        loadComp_MW
        loadComp_MVar
        loadCompIM_powerPercent
        loadCompIM_Rs
        loadCompIM_Rr
        loadCompIM_Lss
        loadCompIM_Lrr
        loadCompIM_Lm
        loadCompIM_R                  % [4 x 4 x numLoadComp]
        loadCompIM_L
        loadCompIM_LTe
        loadCompIM_H
        loadCompIM_m
        loadCompIM_Tm0                % load torque constant, from the init
        loadCompZip_PQabc             % [2 x numLoadComp x 3], (Z,I,P)
        loadCompZip_Tau
        loadCompZip_PQz               % derived at initialisation
        loadCompZip_PQi
        loadCompZip_PQp

        %% GFM IBR
        numIbrGFM = 0
        ibrGFM_bus
        ibrGFM_BaseKV
        ibrGFM_BaseKA
        ibrGFM_MVA_base
        ibrGFM_MW
        ibrGFM_MVar
        ibrGFM_Lf
        ibrGFM_Rf
        ibrGFM_Cf
        ibrGFM_Rcap
        ibrGFM_Lc
        ibrGFM_Rc
        ibrGFM_kC_i
        ibrGFM_kC_p
        ibrGFM_GC
        ibrGFM_kV_i
        ibrGFM_kV_p
        ibrGFM_GV
        ibrGFM_wmeas
        ibrGFM_droopPercentP
        ibrGFM_droopPercentQ
        ibrGFM_Mp                     % derived from the droop percentages
        ibrGFM_Mq
        ibrGFM_w_VR
        ibrGFM_Rv
        ibrGFM_pset                   % set points, from the init
        ibrGFM_qset
        ibrGFM_voDroop_q_set

        %% Triplet stamping
        stampRow
        stampCol
        stampVal
        stampSlice
        map
    end

    methods
        %% ============================ Class ==========================
        function obj = HB(baseMVA,bus,slack_busNum,slack_imScalConst, ...
                numHarmonics,maxDeltaT)
            arguments
                baseMVA
                bus
                slack_busNum
                slack_imScalConst
                numHarmonics
                maxDeltaT
            end

            % Harmonic-balance index globals, used throughout the stamping
            % code. Kept separate from SALT's set because a single session
            % routinely holds both objects.
            hbInitIndices;

            obj.baseMVA = baseMVA;
            obj.bus     = bus;
            obj.numBus  = numel(bus);

            obj.slack_busNum      = slack_busNum;
            obj.slack_imScalConst = slack_imScalConst;

            obj.K = numHarmonics;
            [obj.f0,obj.tSamp,obj.fSamp,obj.N, ...
                obj.GAMMA,obj.invGAMMA,obj.fGamma] = ...
                obj.updateFreqSamples(obj.w0nom);

            obj.maxDeltaT = maxDeltaT;
            [obj.tSim,obj.deltaT] = obj.updateTimeSimSamples;
        end
    end

    methods (Static)
        function hb = fromSalt(salt,numHarmonics,opts)
            %FROMSALT  Build an HB from an already-constructed SALT object.
            %
            %   hb = HB.fromSalt(salt,K) gives a harmonic-balance model of
            %   the SAME system SALT holds, with K harmonics above DC.
            %
            %   WHY NOT A SECOND FILE-DRIVEN CONSTRUCTOR. SALT's constructor
            %   reads the raw case and the device workbook, applies outages
            %   and re-runs the power flow. Duplicating that here would give
            %   two parsers that can drift, and - worse for the comparison
            %   this exists to support - two power-flow solves whose results
            %   could differ in the last digits. Copying SALT's already
            %   parsed parameters makes the two solvers see bit-identical
            %   input data by construction, which is what lets the K = 1
            %   harmonics-off test assert EQUALITY rather than closeness.
            arguments
                salt         (1,1) SALT
                numHarmonics (1,1) double
                opts.maxDeltaT (1,1) double = salt.maxDeltaT
            end

            hb = HB(salt.baseMVA,salt.bus,salt.slack_busNum, ...
                salt.slack_imScalConst,numHarmonics,opts.maxDeltaT);
            hb.caseName   = salt.caseName;
            hb.bus_BaseKv = salt.bus_BaseKv;
            hb.bus_IBase  = salt.bus_IBase;

            % Every network and device PARAMETER SALT parsed. The two
            % classes use the same property names for them - which the
            % composite-load property block already notes is load-bearing -
            % so this is a loop over the shared names rather than 150
            % assignments, and a parameter added to SALT cannot be silently
            % forgotten here.
            %
            % SKIPPED DELIBERATELY: solver state, the harmonic grid, and
            % every OPERATING-POINT quantity. The last group is the subtle
            % one: genrou_efd0, excDc1_vref0, tgov1d_Pref0, loadCompIM_Tm0,
            % the ZIP coefficients and the inverter set points are not
            % parameters - each solver derives them in its own
            % initialisation, and in HB's case some carry the *N of the
            % unnormalized DFT convention, so copying SALT's would be wrong.
            skip = ["baseMVA","bus","numBus","caseName", ...
                    "bus_BaseKv","bus_IBase", ...
                    "slack_busNum","slack_imScalConst", ...
                    "Y","I","x_v","x_v0", ...
                    "maxDeltaT","tSim","deltaT", ...
                    "stampRow","stampCol","stampVal","stampSlice","map", ...
                    "genrou_efd0","genrou_TorqueM0","excDc1_vref0", ...
                    "tgov1d_Pref0","loadCompIM_Tm0", ...
                    "loadCompZip_PQz","loadCompZip_PQi","loadCompZip_PQp", ...
                    "ibrGFM_pset","ibrGFM_qset","ibrGFM_voDroop_q_set"];
            % Constant properties share names across the two classes too
            % (w0, the equation counts, absTolCheck) and cannot be assigned,
            % so the candidate list comes from the metaclass rather than
            % from PROPERTIES.
            pl   = metaclass(hb).PropertyList;
            pHb  = string({pl.Name});
            pHb  = pHb(~[pl.Constant] & strcmp({pl.SetAccess},'public'));
            shared = setdiff(intersect(pHb,string(properties(salt))),skip);
            for p = reshape(shared,1,[])
                hb.(p) = salt.(p);
            end

            % Parameters HB spells differently. The HB name says what the
            % gain IS - Lad/Rfd, the machine-to-exciter field-circuit gain -
            % where SALT's said only which two devices it joins.
            alias = ["excDc1_genToExc","excDc1_genrouToExcDc1"];
            for j = 1:size(alias,1)
                hb.(alias(j,1)) = salt.(alias(j,2));
            end

            % Self-policing. Anything HB can be given that the copy above
            % did not reach has to be here by intent, not by accident: a
            % property renamed on one side only would otherwise be left
            % empty and surface much later as a size mismatch deep in the
            % initialisation. So the leftovers are checked against the list
            % of things HB genuinely owns.
            hbOnly = [ ...
                ... % the harmonic grid, built by the constructor
                "K","N","f0","tSamp","fSamp","fGamma","GAMMA","invGAMMA", ...
                ... % the harmonic source, off unless asked for
                "includeHigherHarmonics","higherHarmonics_bus", ...
                "higherHarmonics_scaleFundamental", ...
                "higherHarmonics_harmonicNum","higherHarmonics_I", ...
                ... % MNA paths SALT does not carry counters for. HB's
                ... % defaults are 0, which is also what BUILDSTAMPMAPS
                ... % insists on, so leaving them is correct.
                "tx_mnaIdxMap","tx_numMnaVars", ...
                "load_mnaIdxMap","load_numMnaVars"];
            leftover = setdiff(setdiff(pHb,shared),[skip,alias(:,1).',hbOnly]);
            assert(isempty(leftover), ...
                ['HB.fromSalt: no source for ' strjoin(leftover,', ') ...
                 '. Either SALT renamed it (add an alias), or it belongs ' ...
                 'in the skip or hbOnly list.'])

            hb = hb.systemValInitializations(salt.pf_busVm,salt.pf_busVa);
            hb = hb.buildStampMaps;
        end
    end

    methods
        %% ==================== Harmonic-balance grid ==================
        function [f0,tSamp,fSamp,N,GAMMA,invGAMMA,fGamma] = ...
                updateFreqSamples(obj,w0)
            %UPDATEFREQSAMPLES  Rebuild the sample grid for a new w0.
            %
            % Called once per Newton iteration, because every tone is
            % k*w0 and w0 is an unknown. Rebuilding here is also what makes
            % the Park transform free of any w0 derivative: the grid is
            % tSamp(n) = n/(N*f0), so w0*tSamp(n) = 2*pi*n/N, independent
            % of w0. The frequency dependence is absorbed into the
            % normalised time variable and survives only in the j*k*w0
            % derivative terms, which ARE differentiated.
            f0 = w0/2/pi;
            frequencies = (0:obj.K).'*f0;

            % All tones are commensurate, so an ordinary DFT suffices.
            [tSamp,fSamp] = hbTimeSamples(obj.K,f0);
            N = numel(tSamp);

            [GAMMA_full,invGAMMA_full] = hbDftMatrices(N);
            fGamma_full = fSamp((floor(N/2)+1):end);
            assert(all(abs(fGamma_full-unique(abs(fSamp)))<obj.absTolCheck))

            % Keep only the retained harmonics. With N = 2*K+1 this
            % selection is the identity; it exists so that oversampling the
            % time grid (see HBTIMESAMPLES) needs no other change.
            [idxReIm,idxFreq] = hbSelectDftHarmonics(fGamma_full,frequencies);
            GAMMA    = GAMMA_full(idxReIm,:);
            invGAMMA = invGAMMA_full(:,idxReIm);
            fGamma   = fGamma_full(idxFreq);
            assert(all(abs(fGamma-frequencies)<obj.absTolCheck))
        end

        function [tSim,deltaT] = updateTimeSimSamples(obj)
            %UPDATETIMESIMSAMPLES  A fine grid, one period, aligned to tSamp.
            delta_tSamp = obj.tSamp(2)-obj.tSamp(1);
            T = obj.tSamp(end) + delta_tSamp;
            frac = ceil(delta_tSamp/obj.maxDeltaT);
            deltaT = delta_tSamp/frac;
            tSim = (0:deltaT:T).';
            assert(tSim(end)==T)
            assert(all(ismembertol(obj.tSamp,tSim)))
        end

        function x = getTimeDomainSignal_iDFT(obj,tSamp,funGetIndices)
            x = hbIdft(tSamp,funGetIndices,obj.x_v,(0:obj.K)*obj.f0,obj.K+1);
        end
    end

    methods
        %% ========================= Index helpers =====================
        % ------------------------- x_v layout ------------------------
        % Two conventions, each chosen for the structure it exposes.
        %
        % NETWORK - harmonic major. Harmonic k occupies a contiguous block
        % of getLinSysSize rows, and within it the layout is SALT's:
        %
        %   [ phase A: bus 1 (re,im) ... bus Nb, mna 1 ... mna Nm
        %     phase B: same
        %     phase C: same ]
        %
        % The network is LTI, so its Jacobian is block diagonal across
        % these harmonic blocks, and keeping them contiguous makes that
        % visible.
        %
        % DEVICES - device major, harmonic inner. One device's COMPLETE
        % harmonic-coupled block is contiguous:
        %
        %   base + k*nEqns + localIndex,   k = 0..K
        %
        % That is what lets SALT's LOCALPATTERN / DEVICETRIPLETS /
        % BLOCKJX / WRITEBLOCKRHS carry over unchanged with nEqns replaced
        % by nEqns*(K+1): the harmonic expansion becomes part of the local
        % pattern rather than something the replication has to know about.
        %
        % Then the GENROU blocks, the EXDC1 blocks, the TGOV1D blocks, and
        % finally one single row for w0.
        %
        % This differs from the legacy EmtHb, which ordered devices
        % type-major then harmonic. validation/validateHb.m carries the
        % permutation between the two.

        function rowIdxAbc = getRowIdxAbc_fromBusNums(obj,busNum,k)
            [~,bus_idx] = ismembertol(busNum,obj.bus);
            assert(size(k,1)==1)
            nNode = obj.numBus + obj.getNumMnaVars;
            idxReIm = obj.numReImEqns*(reshape(bus_idx,1,[])-1) + ...
                (1:obj.numReImEqns).';
            rowIdxAbc = obj.getLinSysSize*k + ...
                [idxReIm(:); ...
                 obj.numReImEqns*nNode   + idxReIm(:); ...
                 2*obj.numReImEqns*nNode + idxReIm(:)];
        end

        function rowIdxAbc = getRowIdxAbc_fromBusNum(obj,busNum,k)
            rowIdxAbc = obj.getRowIdxAbc_fromBusNums(busNum,k);
        end

        function numExtraVars = getNumMnaVars(obj)
            numExtraVars = obj.tx_numMnaVars + obj.xfmr_numMnaVars + ...
                obj.load_numMnaVars;
            assert(numExtraVars == size(obj.tx_mnaIdxMap,1) + ...
                size(obj.xfmr_mnaIdxMap,1) + size(obj.load_mnaIdxMap,1))
        end

        function rowIdxAbc = getRowIdxAbc_MNA_fromMnaVarIdx(obj,varIdx,k)
            assert(size(k,1)==1)
            nNode = obj.numBus + obj.getNumMnaVars;
            idxReIm = obj.numReImEqns*obj.numBus + ...
                obj.numReImEqns*(reshape(varIdx,1,[])-1) + ...
                (1:obj.numReImEqns).';
            rowIdxAbc = obj.getLinSysSize*k + ...
                [idxReIm(:); ...
                 obj.numReImEqns*nNode   + idxReIm(:); ...
                 2*obj.numReImEqns*nNode + idxReIm(:)];
        end

        function obj = addMnaVar(obj,deviceIdx,deviceName)
            numNew = numel(deviceIdx);
            varIdx = obj.getNumMnaVars + (1:numNew).';
            mapName = string(deviceName) + "_mnaIdxMap";
            obj.(mapName) = [obj.(mapName); [deviceIdx(:) varIdx]];
            countName = string(deviceName) + "_numMnaVars";
            obj.(countName) = obj.(countName) + numNew;
        end

        function rowIdx = getRowIdx_fromGenrouIndices(obj,dev,k)
            assert(size(k,1)==1)
            rowIdx = obj.getLinSysSizeAll + ...
                obj.numGenrouEqns*(obj.K+1)*(reshape(dev,[],1)-1) + ...
                obj.numGenrouEqns*k + 1;
        end

        function rowIdx = getRowIdx_fromExcDc1Indices(obj,dev,k)
            assert(size(k,1)==1)
            rowIdx = obj.getLinSysSizeAll + ...
                obj.numGenrouEqns*(obj.K+1)*obj.numGenrou + ...
                obj.numExcDc1Eqns*(obj.K+1)*(reshape(dev,[],1)-1) + ...
                obj.numExcDc1Eqns*k + 1;
        end

        function rowIdx = getRowIdx_fromTgov1dIndices(obj,dev,k)
            assert(size(k,1)==1)
            rowIdx = obj.getLinSysSizeAll + ...
                obj.numGenrouEqns*(obj.K+1)*obj.numGenrou + ...
                obj.numExcDc1Eqns*(obj.K+1)*obj.numExcDc1 + ...
                obj.numTgov1dEqns*(obj.K+1)*(reshape(dev,[],1)-1) + ...
                obj.numTgov1dEqns*k + 1;
        end

        function rowIdx = getRowIdx_fromIbrGFMIndices(obj,dev,k)
            assert(size(k,1)==1)
            rowIdx = obj.getLinSysSizeAll + ...
                obj.numGenrouEqns*(obj.K+1)*obj.numGenrou + ...
                obj.numExcDc1Eqns*(obj.K+1)*obj.numExcDc1 + ...
                obj.numTgov1dEqns*(obj.K+1)*obj.numTgov1d + ...
                obj.numIbrGFMEqns*(obj.K+1)*(reshape(dev,[],1)-1) + ...
                obj.numIbrGFMEqns*k + 1;
        end

        function rowIdx = getRowIdx_fromLoadCompIndices(obj,dev,k)
            assert(size(k,1)==1)
            rowIdx = obj.getLinSysSizeAll + ...
                obj.numGenrouEqns*(obj.K+1)*obj.numGenrou + ...
                obj.numExcDc1Eqns*(obj.K+1)*obj.numExcDc1 + ...
                obj.numTgov1dEqns*(obj.K+1)*obj.numTgov1d + ...
                obj.numIbrGFMEqns*(obj.K+1)*obj.numIbrGFM + ...
                obj.numLoadCompEqns*(obj.K+1)*(reshape(dev,[],1)-1) + ...
                obj.numLoadCompEqns*k + 1;
        end

        function linSysSize = getLinSysSize(obj)
            % Network unknowns at ONE harmonic
            linSysSize = 3*obj.numReImEqns*(obj.numBus + obj.getNumMnaVars);
        end

        function linSysSizeAll = getLinSysSizeAll(obj)
            linSysSizeAll = obj.getLinSysSize*(obj.K+1);
        end

        function sysSizeAll = getSysSizeAll(obj)
            % Device order in x_v: GENROU, EXDC1, TGOV1D, GFM IBR,
            % composite load, then the single w0 row.
            sysSizeAll = obj.getLinSysSizeAll + (obj.K+1)*( ...
                obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + ...
                obj.numTgov1dEqns*obj.numTgov1d + ...
                obj.numIbrGFMEqns*obj.numIbrGFM + ...
                obj.numLoadCompEqns*obj.numLoadComp) + ...
                1;                                      % w0
        end

        function rowIdx = getRowIdx_w0(obj)
            rowIdx = obj.getSysSizeAll;
            if ~isempty(obj.x_v)
                assert(rowIdx==numel(obj.x_v))
            end
        end

        %% ==================== Per-unit conversions ===================
        % Bus is line-to-line RMS kV/kA. The machine base is
        % line-to-neutral RMS, which differs from Kundur (peak); because
        % Vbase_gen = Vbase_bus*1000*sqrt(2/3), the voltage factor between
        % the two works out to exactly 1.
        function c = busToGen_I(obj,bus_idx,gen_idx)
            c = obj.bus_IBase(bus_idx)*1000./obj.genrou_IsBaseRMS(gen_idx);
        end
        function c = genToBus_I(obj,bus_idx,gen_idx)
            c = obj.genrou_IsBaseRMS(gen_idx)./(obj.bus_IBase(bus_idx)*1000);
        end
        function c = busToGen_V(~,~,~), c = 1; end
        function c = genToBus_V(~,~,~), c = 1; end

        % Composite load. Same derivation as the machine, so the voltage
        % factor is again exactly 1; the current needs the RMS motor base.
        function c = busToLoadComp_I(obj,bus_idx,idx)
            c = obj.bus_IBase(bus_idx)*1000./obj.loadComp_IsBaseRMS(idx);
        end
        function c = loadCompToBus_I(obj,bus_idx,idx)
            c = obj.loadComp_IsBaseRMS(idx)./(obj.bus_IBase(bus_idx)*1000);
        end
        function c = busToLoadComp_V(~,~,~), c = 1; end
        function c = loadCompToBus_V(~,~,~), c = 1; end

        % GFM IBR. Unlike the machine, the inverter's per-unit voltage is
        % NOT referred to the bus base, so the voltage factor is a real
        % conversion: Vpu_bus = Vpu_ibr*Vbase_ibr/Vbase_bus.
        function c = busToIbrGFM_I(obj,bus_idx,idx)
            c = obj.bus_IBase(bus_idx)./obj.ibrGFM_BaseKA(idx);
        end
        function c = ibrGFMToBus_I(obj,bus_idx,idx)
            c = obj.ibrGFM_BaseKA(idx)./obj.bus_IBase(bus_idx);
        end
        function c = busToIbrGFM_V(obj,bus_idx,idx)
            c = obj.bus_BaseKv(bus_idx)./obj.ibrGFM_BaseKV(idx);
        end
        function c = ibrGFMToBus_V(obj,bus_idx,idx)
            c = obj.ibrGFM_BaseKV(idx)./obj.bus_BaseKv(bus_idx);
        end
    end

    methods
        %% ====================== Device parameters ====================
        function obj = setTxParams(obj,busFrom,busTo,R,L,C)
            obj.tx_busFrom = busFrom(:);
            obj.tx_busTo   = busTo(:);
            obj.tx_R = R(:); obj.tx_L = L(:); obj.tx_C = C(:);
            % A lossless element has no finite admittance at DC, so it gets
            % an explicit branch-current unknown (MNA) instead.
            idx = find(obj.tx_R==0);
            if ~isempty(idx), obj = obj.addMnaVar(idx,"tx"); end
        end

        function obj = setXfmrParams(obj,busFrom,busTo,R,L)
            obj.xfmr_busFrom = busFrom(:);
            obj.xfmr_busTo   = busTo(:);
            obj.xfmr_R = R(:); obj.xfmr_L = L(:);
            idx = find(obj.xfmr_R==0);
            if ~isempty(idx), obj = obj.addMnaVar(idx,"xfmr"); end
        end

        function obj = setLoadParams(obj,load_bus,load_Z)
            obj.load_bus = load_bus(:);
            obj.load_Z   = load_Z(:);
            obj.load_R   = real(obj.load_Z);
            obj.load_X   = imag(obj.load_Z);
            idx = find((obj.load_R==0) & (obj.load_X>0));
            if ~isempty(idx), obj = obj.addMnaVar(idx,"load"); end
        end

        function obj = setGenrouParams(obj,bus, EsBase,IsBase,BaseMVA, ...
                MW0,Mvar0, Ra,Rfd,R1d,R1q,R2q, Lad,Laq,L0,Ll, ...
                Lffd,Lf1d,L11d,L11q,L22q, H,KD, Xd,Xq)
            obj.genrou_bus = bus(:);
            obj.numGenrou  = numel(bus);

            obj.genrou_EsBase    = EsBase(:);
            obj.genrou_IsBase    = IsBase(:);
            % The bus is RMS and the machine base is peak, so the current
            % conversion needs the RMS machine base.
            obj.genrou_IsBaseRMS = IsBase(:)/sqrt(2);
            obj.genrou_BaseMVA   = BaseMVA(:);
            obj.genrou_MW0       = MW0(:);
            obj.genrou_Mvar0     = Mvar0(:);
            obj.genrou_H         = H(:);
            obj.genrou_KD        = KD(:);

            % Reused from SALT: identical winding order and identical
            % values, verified entry by entry against the legacy
            % setGenParams_matrices.
            [obj.genrou_RMatrix,obj.genrou_LMatrix,obj.genrou_LspMatrix] = ...
                setGenrouParams_matrices(Ra,Rfd,R1d,R1q,R2q, ...
                Lad,Laq,Ll,L0, Lffd,Lf1d,L11d,L11q,L22q);

            obj.genrou_Lad = Lad(:);
            obj.genrou_Rfd = Rfd(:);
            obj.genrou_Xd  = Xd(:);
            obj.genrou_Xq  = Xq(:);
            obj.genrou_Ra  = Ra(:);
        end

        function obj = setExcDc1Params(obj,bus,Tr,Ta,Tc,Tb,Te,Tf1,Kf1,Ka,Ke)
            obj.excDc1_bus = bus(:);
            obj.numExcDc1  = numel(bus);
            if obj.numExcDc1==0, return, end
            obj.excDc1_Tr = Tr(:);  obj.excDc1_Ta = Ta(:);
            obj.excDc1_Tc = Tc(:);  obj.excDc1_Tb = Tb(:);
            obj.excDc1_Te = Te(:);  obj.excDc1_Tf1 = Tf1(:);
            obj.excDc1_Kf1 = Kf1(:); obj.excDc1_Ka = Ka(:);
            obj.excDc1_Ke  = Ke(:);
            [found,gi] = ismembertol(obj.excDc1_bus,obj.genrou_bus);
            assert(all(found),'Every EXDC1 must sit on a GENROU bus')
            obj.excDc1_genrouIdx = gi(:);
            % Field-circuit gain: efd = (Lad/Rfd)*x1
            obj.excDc1_genToExc = obj.genrou_Lad(gi)./obj.genrou_Rfd(gi);
        end

        function obj = setTgov1dParams(obj,bus,T1,T2,T3,R,Dt)
            obj.tgov1d_bus = bus(:);
            obj.numTgov1d  = numel(bus);
            if obj.numTgov1d==0, return, end
            obj.tgov1d_T1 = T1(:); obj.tgov1d_T2 = T2(:);
            obj.tgov1d_T3 = T3(:); obj.tgov1d_R  = R(:);
            obj.tgov1d_Dt = Dt(:);
            [found,gi] = ismembertol(obj.tgov1d_bus,obj.genrou_bus);
            assert(all(found),'Every TGOV1D must sit on a GENROU bus')
            obj.tgov1d_genrouIdx = gi(:);
        end

        % ***** COMPOSITE LOAD *****
        % Copied from SALT.setLoadCompParams unchanged: it is pure
        % parameter assembly with no single-harmonic assumption in it.
        function obj = setLoadCompParams(obj,bus, ...
                VBase,IBase,MVA_base, MW,MVar, ...
                IM_powerPercent, IM_Rs,IM_Rr, IM_Lss,IM_Lrr,IM_Lm, ...
                IM_H, IM_m, Zip_PQabc, Zip_Tau)

            obj.loadComp_bus = bus(:);
            obj.numLoadComp  = numel(obj.loadComp_bus);
            obj.loadComp_VBase = VBase(:);
            obj.loadComp_IBase = IBase(:);
            obj.loadComp_IsBaseRMS = IBase(:)/sqrt(2);
            obj.loadComp_MVA_base = MVA_base(:);
            obj.loadComp_MW = MW(:);
            obj.loadComp_MVar = MVar(:);
            obj.loadCompIM_powerPercent = IM_powerPercent(:);
            obj.loadCompIM_Rs = IM_Rs(:);   obj.loadCompIM_Rr = IM_Rr(:);
            obj.loadCompIM_Lss = IM_Lss(:); obj.loadCompIM_Lrr = IM_Lrr(:);
            obj.loadCompIM_Lm = IM_Lm(:);

            n = obj.numLoadComp;
            % (R) stator and rotor resistance, block diagonal in dq
            obj.loadCompIM_R = zeros(4,4,n);
            obj.loadCompIM_R(1,1,:) = IM_Rs;
            obj.loadCompIM_R(2,2,:) = IM_Rs;
            obj.loadCompIM_R(3,3,:) = IM_Rr;
            obj.loadCompIM_R(4,4,:) = IM_Rr;
            % (L) self and mutual inductance
            obj.loadCompIM_L = zeros(4,4,n);
            obj.loadCompIM_L(1,1,:) = IM_Lss;
            obj.loadCompIM_L(2,2,:) = IM_Lss;
            obj.loadCompIM_L(3,3,:) = IM_Lrr;
            obj.loadCompIM_L(4,4,:) = IM_Lrr;
            obj.loadCompIM_L(1,3,:) = IM_Lm;
            obj.loadCompIM_L(2,4,:) = IM_Lm;
            obj.loadCompIM_L(3,1,:) = IM_Lm;
            obj.loadCompIM_L(4,2,:) = IM_Lm;
            % (Swing) electrical torque operator
            matSelect = [0 0; 0 0; 0 -1; 1 0]*[0 0 1 0;0 0 0 1];
            obj.loadCompIM_LTe = pagemtimes( ...
                pagetranspose(obj.loadCompIM_L),matSelect);
            obj.loadCompIM_m = IM_m(:);
            obj.loadCompIM_H = IM_H(:);

            % ZIP shares must sum to one within P and within Q
            assert(isequal(size(Zip_PQabc),[2,n,3]))
            assert(all(abs(sum(Zip_PQabc,3)-1)<obj.absTolCheck,'all'))
            obj.loadCompZip_PQabc = Zip_PQabc;
            assert(all(Zip_Tau>0))
            obj.loadCompZip_Tau = Zip_Tau;
        end

        % ***** GFM IBR *****
        % Copied from SALT.setIbrGFMParams unchanged, plus the droop gains,
        % which SALT derives in its own initialisation.
        function obj = setIbrGFMParams(obj,bus, ...
                BaseKV,BaseKA,MVA_base, MW,MVar, ...
                Lf,Rf,Cf,Rcap,Lc,Rc, kC_i,kC_p,GC, kV_i,kV_p,GV, ...
                wmeas,droopPercentP,droopPercentQ,w_VR,Rv)

            obj.ibrGFM_bus = bus(:);
            obj.numIbrGFM  = numel(obj.ibrGFM_bus);

            obj.ibrGFM_BaseKV = BaseKV(:);
            obj.ibrGFM_BaseKA = BaseKA(:);
            obj.ibrGFM_MVA_base = MVA_base(:);
            % Check the sqrt(3), since the base assumes V line-to-line, RMS
            assert(all(abs(BaseKV(:).*BaseKA(:)*sqrt(3)-MVA_base(:)) < ...
                obj.absTolCheck))

            obj.ibrGFM_MW = MW(:);  obj.ibrGFM_MVar = MVar(:);
            obj.ibrGFM_Lf = Lf(:); obj.ibrGFM_Rf = Rf(:);
            obj.ibrGFM_Cf = Cf(:); obj.ibrGFM_Rcap = Rcap(:);
            obj.ibrGFM_Lc = Lc(:); obj.ibrGFM_Rc = Rc(:);
            obj.ibrGFM_kC_i = kC_i(:); obj.ibrGFM_kC_p = kC_p(:);
            obj.ibrGFM_GC = GC(:);
            obj.ibrGFM_kV_i = kV_i(:); obj.ibrGFM_kV_p = kV_p(:);
            obj.ibrGFM_GV = GV(:);
            obj.ibrGFM_wmeas = wmeas(:);
            obj.ibrGFM_droopPercentP = droopPercentP(:);
            obj.ibrGFM_droopPercentQ = droopPercentQ(:);
            obj.ibrGFM_w_VR = w_VR(:); obj.ibrGFM_Rv = Rv(:);
            % Droop gains, exactly as SALT derives them in its own
            % initialisation (Mp = droopPercentP*w0, Mq = droopPercentQ).
            % Done here instead because they depend only on parameters, not
            % on any solved state.
            obj.ibrGFM_Mp = droopPercentP(:)*obj.w0nom;
            obj.ibrGFM_Mq = droopPercentQ(:);
        end

        function obj = setHigherHarmonicsParams(obj,bus,scaleFundamental, ...
                harmonicNum,Ipf)
            %SETHIGHERHARMONICSPARAMS  A harmonic current source at a bus.
            %
            % Ipf is the fundamental current at that bus from the power
            % flow; the injected source is scaleFundamental times it,
            % placed at harmonicNum and spread over three balanced phases.
            %
            % NOTE: scaleFundamental is a BARE MULTIPLIER. The legacy name
            % said "percent" but the arithmetic never divided by 100, so a
            % value of 2 is 200% of the fundamental, not 2%.
            obj.includeHigherHarmonics = true;
            obj.higherHarmonics_bus = bus(:);
            obj.higherHarmonics_scaleFundamental = scaleFundamental(:);
            obj.higherHarmonics_harmonicNum = harmonicNum(:);

            nH = numel(bus);
            obj.higherHarmonics_I = zeros(3*obj.numReImEqns,nH);
            % x_v holds unnormalized coefficients (see HBIDFT), hence the
            % *N; the /2 is the split onto the conjugate pair, which only
            % applies away from DC.
            assert(all(harmonicNum~=0), ...
                'A DC harmonic source would not take the 1/2 conjugate split')
            Ipf_dft = Ipf(:)*obj.N/2;
            for j = 1:nH
                [IphA,IphB,IphC] = calc3Ph( ...
                    obj.higherHarmonics_scaleFundamental(j)*Ipf_dft(j));
                obj.higherHarmonics_I(:,j) = ...
                    [real(IphA);imag(IphA); ...
                     real(IphB);imag(IphB); ...
                     real(IphC);imag(IphC)];
            end
        end
    end

    methods
        %% ======================= Initialization ======================
        function obj = systemValInitializations(obj,bus_Vm,bus_Va)
            %SYSTEMVALINITIALIZATIONS  Seed x_v from the power flow.
            %
            % The entire initial guess is the power-flow operating point
            % placed at the FUNDAMENTAL only - bus voltages balanced
            % three-phase at k = 1, machine states at k = 0, every other
            % harmonic zero. So the starting point is "the single-harmonic
            % answer with no harmonic content", and Newton has to grow the
            % harmonics from the injected source.
            global idxHbGenrou_Idq0 idxHbGenrou_Edq0 idxHbGenrou_Efd ...
                   idxHbGenrou_Wr idxHbGenrou_DeltaTheta ...
                   idxHbGenrou_IphABC idxHbGenrou_Tm                %#ok<GVMIS>
            global idxHbExcDc1_x1 idxHbExcDc1_x2 idxHbExcDc1_x3 ...
                   idxHbExcDc1_x4 idxHbExcDc1_x5 idxHbExcDc1_Vt    %#ok<GVMIS>
            global idxHbTgov1d_x1 idxHbTgov1d_x2 idxHbTgov1d_Pmech %#ok<GVMIS>
            global idxHbIbrGFM_IphABC idxHbIbrGFM_Vg idxHbIbrGFM_Vo ...
                   idxHbIbrGFM_If idxHbIbrGFM_Io idxHbIbrGFM_Vi ...
                   idxHbIbrGFM_Gamma idxHbIbrGFM_IfRef idxHbIbrGFM_X ...
                   idxHbIbrGFM_Wr idxHbIbrGFM_VoRefDroop_q ...
                   idxHbIbrGFM_PqAvg idxHbIbrGFM_DeltaTheta ...
                   idxHbIbrGFM_Vo_VR_avg idxHbIbrGFM_Vo_VR ...
                   idxHbIbrGFM_VoRef idxHbIbrGFM_Io0                %#ok<GVMIS>
            global idxHbLoadComp_Vs idxHbLoadComp_IphABC ...
                   idxHbLoadComp_WSlip idxHbLoadComp_Is ...
                   idxHbLoadComp_Ir idxHbLoadComp_Wr ...
                   idxHbLoadComp_PQ idxHbLoadComp_Vm ...
                   idxHbLoadComp_Izip                               %#ok<GVMIS>

            obj.x_v = zeros(obj.getSysSizeAll,1);

            %% ---- Linear: bus voltages at the fundamental -------------
            % Reused from SALT; bus voltage is line-to-line, RMS.
            [bus_V3ph,bus_VphA] = busValInitializations(bus_Vm,bus_Va);
            idxBus = obj.getRowIdxAbc_fromBusNums(obj.bus,1);
            obj.x_v(idxBus(1:2:end)) = real(bus_V3ph);
            obj.x_v(idxBus(2:2:end)) = imag(bus_V3ph);

            %% ---- MNA branch currents at the fundamental --------------
            % Only the lossless-transformer path is exercised; the
            % transmission-line and load MNA paths were never debugged in
            % the legacy and are rejected outright in BUILDSTAMPMAPS
            % rather than silently producing a wrong answer.
            wk = obj.w0nom;
            for j = 1:obj.xfmr_numMnaVars
                dev = obj.xfmr_mnaIdxMap(j,1);
                mna = obj.xfmr_mnaIdxMap(j,2);
                idxF = obj.getRowIdxAbc_fromBusNum(obj.xfmr_busFrom(dev),1);
                idxT = obj.getRowIdxAbc_fromBusNum(obj.xfmr_busTo(dev),1);
                VF = complex(obj.x_v(idxF(1:2:end)),obj.x_v(idxF(2:2:end)));
                VT = complex(obj.x_v(idxT(1:2:end)),obj.x_v(idxT(2:2:end)));
                assert(obj.xfmr_R(dev)==0)
                I_FT = (VF-VT)*HB.calcYeq_RLSeries(0,obj.xfmr_L(dev),wk);
                idxI = obj.getRowIdxAbc_MNA_fromMnaVarIdx(mna,1);
                obj.x_v(idxI(1:2:end)) = real(I_FT);
                obj.x_v(idxI(2:2:end)) = imag(I_FT);
            end

            %% ---- Nonlinear devices, at DC ----------------------------
            genrou_VAbsPhA = []; genrou_efd = []; genrou_TorqueM = [];
            if obj.numGenrou > 0
                [id,iq,ifd,ed,eq,efd,THETA,TorqueM,genrou_VAbsPhA, ...
                    IComplexPhA_busPU] = ...
                    hbGenrouValInitCalculations(obj,bus_VphA,bus_Va);
                genrou_efd = efd; genrou_TorqueM = TorqueM;
                r = obj.getRowIdx_fromGenrouIndices((1:obj.numGenrou).',0)-1;
                % Only real parts: these all sit at DC.
                obj.x_v(r+idxHbGenrou_Idq0(1))  = id;
                obj.x_v(r+idxHbGenrou_Idq0(3))  = iq;
                obj.x_v(r+idxHbGenrou_Idq0(7))  = ifd;
                obj.x_v(r+idxHbGenrou_Edq0(1))  = ed;
                obj.x_v(r+idxHbGenrou_Edq0(3))  = eq;
                obj.x_v(r+idxHbGenrou_Efd(1))   = efd;
                obj.x_v(r+idxHbGenrou_Wr(1))    = obj.w0nom;
                obj.x_v(r+idxHbGenrou_DeltaTheta(1)) = THETA;
                obj.x_v(r+idxHbGenrou_Tm(1))    = TorqueM;
                % Terminal current, at the FUNDAMENTAL. The dq states above
                % are DC quantities in the rotor frame, but the abc
                % terminal current is a 60 Hz quantity, so it is seeded at
                % k = 1 as the balanced three-phase set implied by the
                % power-flow current, in BUS per-unit. (Seeding it at DC
                % instead - say from the inverse Park transform at t = 0 -
                % would put a real 60 Hz current in the DC bin.)
                r1 = obj.getRowIdx_fromGenrouIndices( ...
                    (1:obj.numGenrou).',1)-1;
                for g = 1:obj.numGenrou
                    [IphA,IphB,IphC] = calc3Ph(IComplexPhA_busPU(g));
                    obj.x_v(r1(g)+idxHbGenrou_IphABC) = [ ...
                        real(IphA);imag(IphA); ...
                        real(IphB);imag(IphB); ...
                        real(IphC);imag(IphC)];
                end
                % The initialization works in the time domain, so convert
                % the control set points to the unnormalized DFT
                % convention. At DC there is no factor of 2.
                obj.genrou_efd0    = efd * obj.N;
                obj.genrou_TorqueM0 = TorqueM * obj.N;
            end

            if obj.numExcDc1 > 0
                et0 = genrou_VAbsPhA(obj.excDc1_genrouIdx);
                [x1,x2,x3,x4,x5,vref] = hbExcDc1ValInitCalculations( ...
                    obj,et0,genrou_efd(obj.excDc1_genrouIdx));
                r = obj.getRowIdx_fromExcDc1Indices((1:obj.numExcDc1).',0)-1;
                obj.x_v(r+idxHbExcDc1_x1(1)) = x1;
                obj.x_v(r+idxHbExcDc1_x2(1)) = x2;
                obj.x_v(r+idxHbExcDc1_x3(1)) = x3;
                obj.x_v(r+idxHbExcDc1_x4(1)) = x4;
                obj.x_v(r+idxHbExcDc1_x5(1)) = x5;
                obj.x_v(r+idxHbExcDc1_Vt(1)) = et0;
                obj.excDc1_vref0 = vref * obj.N;
            end

            if obj.numTgov1d > 0
                [t1,t2,pmech,pref] = hbTgov1dValInitCalculations( ...
                    obj,genrou_TorqueM(obj.tgov1d_genrouIdx));
                r = obj.getRowIdx_fromTgov1dIndices((1:obj.numTgov1d).',0)-1;
                obj.x_v(r+idxHbTgov1d_x1(1))    = t1;
                obj.x_v(r+idxHbTgov1d_x2(1))    = t2;
                obj.x_v(r+idxHbTgov1d_Pmech(1)) = pmech;
                obj.tgov1d_Pref0 = pref * obj.N;
            end

            % A dq pair occupies four slots, (d_re,d_im,q_re,q_im), so the
            % two REAL slots of a pair are its 1st and 3rd. Every value
            % placed below is a DC quantity, hence real.
            dqRe = @(loc) loc([1 3]);

            if obj.numIbrGFM > 0
                % SALT's own routine, used verbatim. It applies the GENUINE
                % 2x3 PSCAD Park to the three-phase phasors at t = 0, so
                % every dq pair it returns is already in the frame the HB
                % inverter kernels use and lands at DC unchanged - no
                % conjugation or sign fix is needed anywhere here.
                [vg_dq,vo_dq,vi_dq,if_dq,io_dq,gamma_dq,ifRef_dq,x_dq, ...
                    wr_ibr,voRef_dq,pqavg,iabc_reIm, ...
                    pset,qset,voDroopQset,THETA_ibr, ...
                    voVRavg_dq,voVR_dq,voDroopQref] = ...
                    ibrGFMValInitCalculations_assumeW0(obj,bus_VphA,bus_Va);
                % Droop set points stay as PHYSICAL values: unlike efd0,
                % vref0 and Pref0, the inverter kernel applies the *N of the
                % unnormalized DFT convention itself when it builds the
                % constant term.
                obj.ibrGFM_pset = pset;
                obj.ibrGFM_qset = qset;
                obj.ibrGFM_voDroop_q_set = voDroopQset;

                r = obj.getRowIdx_fromIbrGFMIndices( ...
                    (1:obj.numIbrGFM).',0)-1;
                obj.x_v(r+dqRe(idxHbIbrGFM_Vg))    = vg_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_Vo))    = vo_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_Vi))    = vi_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_If))    = if_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_Io))    = io_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_Gamma)) = gamma_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_IfRef)) = ifRef_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_X))     = x_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_PqAvg)) = pqavg.';
                obj.x_v(r+dqRe(idxHbIbrGFM_Vo_VR_avg)) = voVRavg_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_Vo_VR))     = voVR_dq.';
                obj.x_v(r+dqRe(idxHbIbrGFM_VoRef))     = voRef_dq.';
                obj.x_v(r+idxHbIbrGFM_Wr(1)) = wr_ibr(:);
                obj.x_v(r+idxHbIbrGFM_VoRefDroop_q(1)) = voDroopQref(:);
                % DeltaTheta, not Theta: HB writes theta = w0*t + deltaTheta
                % and SALT's THETA is theta at t = 0, so the deviation
                % starts at THETA, exactly as the GENROU's does.
                obj.x_v(r+idxHbIbrGFM_DeltaTheta(1)) = THETA_ibr(:);
                % Terminal current at the FUNDAMENTAL, for the same reason
                % as the GENROU's: the dq states are DC in the inverter's
                % own frame, but the abc current is a 60 Hz quantity. SALT
                % already returns it as interleaved (re,im) in BUS per-unit,
                % which is this block's IphABC layout.
                r1 = obj.getRowIdx_fromIbrGFMIndices( ...
                    (1:obj.numIbrGFM).',1)-1;
                obj.x_v(r1+idxHbIbrGFM_IphABC) = iabc_reIm.';
            end

            if obj.numLoadComp > 0
                % SALT's own routine, used verbatim. Its dq pairs are the
                % phase-A PHASOR written as (real,imag) with theta_sys = 0,
                % which is exactly what the composite load's -sin Park gives
                % at DC - see the Park-convention note in
                % CALCHBLOADCOMPMATRICES_FD for why that convention and not
                % the inverter's.
                [vdq_s,loadComp_IphABC,wslip,thetaSys,wr_im, ...
                    idq_s,idq_r,PQ_Zip,Vm_zip,IzipPhA, ...
                    T0,PQz,PQi,PQp] = ...
                    compLoadValInitCalculations_assumeW0(obj,bus_VphA);
                % Load torque constant and ZIP polynomial coefficients.
                % Both are PHYSICAL constants: the kernels evaluate the
                % terms they appear in pointwise in TIME and then transform,
                % so the DFT supplies the N scaling on its own.
                obj.loadCompIM_Tm0   = T0;
                obj.loadCompZip_PQz  = PQz;
                obj.loadCompZip_PQi  = PQi;
                obj.loadCompZip_PQp  = PQp;
                assert(all(abs(thetaSys)<obj.absTolCheck), ...
                    'HB: the composite load frame angle must start at zero')

                r = obj.getRowIdx_fromLoadCompIndices( ...
                    (1:obj.numLoadComp).',0)-1;
                obj.x_v(r+dqRe(idxHbLoadComp_Vs)) = vdq_s.';
                obj.x_v(r+dqRe(idxHbLoadComp_Is)) = idq_s.';
                obj.x_v(r+dqRe(idxHbLoadComp_Ir)) = idq_r.';
                obj.x_v(r+dqRe(idxHbLoadComp_PQ)) = PQ_Zip.';
                obj.x_v(r+dqRe(idxHbLoadComp_Izip)) = ...
                    [real(IzipPhA(:)),imag(IzipPhA(:))];
                obj.x_v(r+idxHbLoadComp_WSlip(1)) = wslip(:);
                obj.x_v(r+idxHbLoadComp_Wr(1))    = wr_im(:);
                obj.x_v(r+idxHbLoadComp_Vm(1))    = Vm_zip(:);
                % DeltaThetaSys stays at zero: the frame rotates at exactly
                % the fundamental, so the deviation is zero at every
                % harmonic, which is SALT's thetaSys = 0 asserted above.
                r1 = obj.getRowIdx_fromLoadCompIndices( ...
                    (1:obj.numLoadComp).',1)-1;
                iphReIm = zeros(6,obj.numLoadComp);
                iphReIm(1:2:end,:) = real(loadComp_IphABC);
                iphReIm(2:2:end,:) = imag(loadComp_IphABC);
                obj.x_v(r1+idxHbLoadComp_IphABC) = iphReIm.';
            end

            %% ---- DFT normalization -----------------------------------
            % Everything above was written as a TIME-DOMAIN value, but x_v
            % holds unnormalized DFT coefficients under
            %     x(n) = (1/N)*sum_k X(k)*exp(1j*2*pi*k*n/N).
            % Two steps follow from that single convention:
            %
            % 1. multiply by N, to undo the 1/N in the inverse transform
            obj.x_v = obj.x_v * obj.N;
            % 2. halve every non-DC component, because a real signal splits
            %    it across a conjugate pair. DC has no partner, so put it
            %    back.
            obj.x_v = obj.x_v / 2;
            obj.x_v(obj.dcRowIdx) = obj.x_v(obj.dcRowIdx) * 2;

            %% ---- Frequency -------------------------------------------
            obj.x_v(obj.getRowIdx_w0) = obj.w0nom;

            obj.x_v0 = obj.x_v;
        end

        function idx = dcRowIdx(obj)
            %DCROWIDX  Every row of x_v that holds a k = 0 coefficient.
            idx = [obj.getRowIdxAbc_fromBusNums(obj.bus,0); ...
                   obj.getRowIdxAbc_MNA_fromMnaVarIdx( ...
                       (1:obj.getNumMnaVars).',0)];
            if obj.numGenrou > 0
                r = obj.getRowIdx_fromGenrouIndices((1:obj.numGenrou).',0)-1;
                idx = [idx; reshape(r.' + (1:obj.numGenrouEqns).',[],1)];
            end
            if obj.numExcDc1 > 0
                r = obj.getRowIdx_fromExcDc1Indices((1:obj.numExcDc1).',0)-1;
                idx = [idx; reshape(r.' + (1:obj.numExcDc1Eqns).',[],1)];
            end
            if obj.numTgov1d > 0
                r = obj.getRowIdx_fromTgov1dIndices((1:obj.numTgov1d).',0)-1;
                idx = [idx; reshape(r.' + (1:obj.numTgov1dEqns).',[],1)];
            end
            if obj.numIbrGFM > 0
                r = obj.getRowIdx_fromIbrGFMIndices((1:obj.numIbrGFM).',0)-1;
                idx = [idx; reshape(r.' + (1:obj.numIbrGFMEqns).',[],1)];
            end
            if obj.numLoadComp > 0
                r = obj.getRowIdx_fromLoadCompIndices((1:obj.numLoadComp).',0)-1;
                idx = [idx; reshape(r.' + (1:obj.numLoadCompEqns).',[],1)];
            end
        end
    end

    methods
        %% ======================= Stamp maps ==========================
        function obj = buildStampMaps(obj)
            %BUILDSTAMPMAPS  Precompute the row/col of every triplet.
            %
            % Nothing here depends on x_v, so it runs once. The stamping
            % methods then only refill obj.stampVal and ASSEMBLE_Y_I turns
            % the three arrays into the sparse Jacobian. Same contract as
            % SALT.buildStampMaps; the harmonic dimension enters in two
            % different ways:
            %
            %   harmonic-DIAGONAL patterns for everything linear - the
            %   network is LTI and the per-harmonic device stamps are
            %   evaluated one k at a time, so harmonic k couples only to k.
            %
            %   harmonic-DENSE patterns for every nonlinear term, because
            %   GAMMA*diag(.)*invGAMMA couples every harmonic to every
            %   other. That mixing is the point of the method.
            global idxHbGenrou_Idq0 idxHbGenrou_Edq0 idxHbGenrou_Efd ...
                   idxHbGenrou_Wr idxHbGenrou_DeltaTheta ...
                   idxHbGenrou_IphABC idxHbGenrou_Tm                %#ok<GVMIS>
            global idxHbExcDc1_x1 idxHbExcDc1_x2 idxHbExcDc1_x3 ...
                   idxHbExcDc1_x4 idxHbExcDc1_x5 idxHbExcDc1_Vt     %#ok<GVMIS>
            global idxHbTgov1d_x1 idxHbTgov1d_x2 idxHbTgov1d_Pmech  %#ok<GVMIS>
            global idxHbIbrGFM_IphABC idxHbIbrGFM_Vg idxHbIbrGFM_Vo ...
                   idxHbIbrGFM_If idxHbIbrGFM_Io idxHbIbrGFM_Vi ...
                   idxHbIbrGFM_Gamma idxHbIbrGFM_IfRef idxHbIbrGFM_X ...
                   idxHbIbrGFM_Wr idxHbIbrGFM_VoRefDroop_q ...
                   idxHbIbrGFM_PqAvg idxHbIbrGFM_DeltaTheta ...
                   idxHbIbrGFM_Vo_VR_avg idxHbIbrGFM_Vo_VR ...
                   idxHbIbrGFM_VoRef idxHbIbrGFM_Io0                %#ok<GVMIS>
            global idxHbLoadComp_Vs idxHbLoadComp_IphABC ...
                   idxHbLoadComp_WSlip idxHbLoadComp_DeltaThetaSys ...
                   idxHbLoadComp_Is idxHbLoadComp_Ir idxHbLoadComp_Wr ...
                   idxHbLoadComp_PQ idxHbLoadComp_Vm ...
                   idxHbLoadComp_Izip                                %#ok<GVMIS>

            % The legacy never debugged the transmission-line and load MNA
            % paths and its frequency-as-a-variable stamp for them errors
            % out, so refuse them here rather than produce a wrong answer.
            assert(obj.tx_numMnaVars==0, ...
                ['HB: a lossless transmission line needs the MNA path, ' ...
                 'which was never completed for frequency as a variable.'])
            assert(obj.load_numMnaVars==0, ...
                'HB: the constant-impedance MNA path is not implemented.')

            % Returns [tSim,deltaT], not obj - the sample grid is a
            % property pair, not a rebuilt object.
            [obj.tSim,obj.deltaT] = obj.updateTimeSimSamples;

            idxW0 = obj.getRowIdx_w0;
            K1 = obj.K+1;
            nG = obj.numGenrou; nE = obj.numExcDc1; nV = obj.numTgov1d;
            nI = obj.numIbrGFM; nL = obj.numLoadComp;
            eqG = obj.numGenrouEqns; eqE = obj.numExcDc1Eqns;
            eqV = obj.numTgov1dEqns; eqI = obj.numIbrGFMEqns;
            eqL = obj.numLoadCompEqns;

            R = {}; C = {}; slice = struct(); n = 0;
            function addSlice(name,r,c)
                R{end+1} = r(:); C{end+1} = c(:);
                slice.(name) = n + (1:numel(r)).';
                n = n + numel(r);
            end

            %% ---------------- Linear network, per harmonic ------------
            % Ordering everywhere below is harmonic-outer, device-inner,
            % entry-innermost, which is the order the stamping loop fills.
            nTx = numel(obj.tx_busFrom);
            rTx = zeros(12*13*nTx*K1,1); cTx = rTx; p = 0;
            rowsTx = zeros(12*nTx,K1);
            for k = 0:obj.K
                for b = 1:nTx
                    idxFT = obj.branchRowIdxFT(obj.tx_busFrom(b), ...
                                               obj.tx_busTo(b),k);
                    [rr,cc] = saltStampTriplets(idxFT,[idxFT;idxW0]);
                    rTx(p+(1:156)) = rr; cTx(p+(1:156)) = cc; p = p+156;
                    rowsTx((b-1)*12+(1:12),k+1) = idxFT;
                end
            end
            addSlice('txY',rTx,cTx);
            obj.map.txIrows = rowsTx;

            % Transformers. A lossless one carries no admittance, so it
            % gets an explicit branch-current unknown instead.
            isMna = false(numel(obj.xfmr_busFrom),1);
            isMna(obj.xfmr_mnaIdxMap(:,1)) = true;
            obj.map.xfmr.isMna  = isMna;
            obj.map.xfmr.idxAdm = find(~isMna);
            obj.map.xfmr.idxMna = find(isMna);

            rY = []; cY = []; rowsXfY = [];
            for k = 0:obj.K
                rowsK = [];
                for b = obj.map.xfmr.idxAdm.'
                    idxFT = obj.branchRowIdxFT(obj.xfmr_busFrom(b), ...
                                               obj.xfmr_busTo(b),k);
                    [rr,cc] = saltStampTriplets(idxFT,[idxFT;idxW0]);
                    rY = [rY;rr]; cY = [cY;cc];                  %#ok<AGROW>
                    rowsK = [rowsK;idxFT];                       %#ok<AGROW>
                end
                rowsXfY = [rowsXfY,rowsK];                       %#ok<AGROW>
            end
            addSlice('xfmrY',rY,cY);
            obj.map.xfmrYIrows = rowsXfY;

            rSrc = []; cSrc = []; rZ = []; cZ = []; rowsXfZ = [];
            for k = 0:obj.K
                rowsK = [];
                for b = obj.map.xfmr.idxMna.'
                    idxFT = obj.branchRowIdxFT(obj.xfmr_busFrom(b), ...
                                               obj.xfmr_busTo(b),k);
                    mna = obj.xfmr_mnaIdxMap(obj.xfmr_mnaIdxMap(:,1)==b,2);
                    idxI = obj.getRowIdxAbc_MNA_fromMnaVarIdx(mna,k);
                    % The branch current is a source in the two node
                    % equations it connects...
                    [rr,cc] = saltStampTriplets(idxFT,idxI);
                    rSrc = [rSrc;rr]; cSrc = [cSrc;cc];          %#ok<AGROW>
                    % ...and gets its own row, V_from - V_to = Z*I
                    [rr,cc] = saltStampTriplets(idxI,[idxFT;idxI;idxW0]);
                    rZ = [rZ;rr]; cZ = [cZ;cc];                  %#ok<AGROW>
                    rowsK = [rowsK;idxI];                        %#ok<AGROW>
                end
                rowsXfZ = [rowsXfZ,rowsK];                       %#ok<AGROW>
            end
            addSlice('xfmrIsrc',rSrc,cSrc);
            addSlice('xfmrZ',rZ,cZ);
            obj.map.xfmrZIrows = rowsXfZ;

            rL = []; cL = []; rowsLd = [];
            for k = 0:obj.K
                rowsK = [];
                for b = 1:numel(obj.load_bus)
                    idxF = obj.getRowIdxAbc_fromBusNum(obj.load_bus(b),k);
                    [rr,cc] = saltStampTriplets(idxF,[idxF;idxW0]);
                    rL = [rL;rr]; cL = [cL;cc];                  %#ok<AGROW>
                    rowsK = [rowsK;idxF];                        %#ok<AGROW>
                end
                rowsLd = [rowsLd,rowsK];                         %#ok<AGROW>
            end
            addSlice('loadY',rL,cL);
            obj.map.loadIrows = rowsLd;

            %% ---------------- Device terminal current -----------------
            % Machines produce power, so their current leaves the device
            % and enters the circuit: stamped negative. Harmonic-diagonal,
            % because it is a plain identity at each k.
            rI = []; cI = [];
            for k = 0:obj.K
                for g = 1:nG
                    vabc = obj.getRowIdxAbc_fromBusNum(obj.genrou_bus(g),k);
                    iabc = obj.getRowIdx_fromGenrouIndices(g,k)-1 + ...
                        idxHbGenrou_IphABC;
                    [rr,cc] = saltStampTriplets(vabc,iabc);
                    rI = [rI;rr]; cI = [cI;cc];                  %#ok<AGROW>
                end
            end
            addSlice('injGenrou',rI,cI);

            % The inverter also produces, so it is stamped negative too;
            % the composite load consumes, so its current flows INTO the
            % device and is stamped positive.
            rI = []; cI = [];
            for k = 0:obj.K
                for g = 1:nI
                    vabc = obj.getRowIdxAbc_fromBusNum(obj.ibrGFM_bus(g),k);
                    iabc = obj.getRowIdx_fromIbrGFMIndices(g,k)-1 + ...
                        idxHbIbrGFM_IphABC;
                    [rr,cc] = saltStampTriplets(vabc,iabc);
                    rI = [rI;rr]; cI = [cI;cc];                  %#ok<AGROW>
                end
            end
            addSlice('injIbrGFM',rI,cI);

            rI = []; cI = [];
            for k = 0:obj.K
                for g = 1:nL
                    vabc = obj.getRowIdxAbc_fromBusNum(obj.loadComp_bus(g),k);
                    iabc = obj.getRowIdx_fromLoadCompIndices(g,k)-1 + ...
                        idxHbLoadComp_IphABC;
                    [rr,cc] = saltStampTriplets(vabc,iabc);
                    rI = [rI;rr]; cI = [cI;cc];                  %#ok<AGROW>
                end
            end
            addSlice('injLoadComp',rI,cI);

            %% ---------------- GENROU ----------------------------------
            gi = idxHbGenrou_Idq0; ge = idxHbGenrou_Edq0;
            gE = idxHbGenrou_Efd;  gw = idxHbGenrou_Wr;
            gt = idxHbGenrou_DeltaTheta; gia = idxHbGenrou_IphABC;
            gT = idxHbGenrou_Tm;

            % Linear: one block list, evaluated once per harmonic. One
            % COLUMN GROUP per entry, so a kernel can hand back each block
            % flattened with no interleaving.
            Plin = { gi(1:8),  [ge gE]       % stator voltage (d,q,0,fd)
                     gi,       gi            % -R + j*wk*L
                     gw,       gw            % swing inertia and damping
                     gw,       gT            % mechanical torque
                     gt,       gw            % d(deltaTheta)/dt = wr - w0
                     gt,       gt            % -j*wk*deltaTheta
                     gia,       gia };         % terminal current scaling
            obj.map.genrouLin = obj.localPatternHbDiag(Plin,eqG);
            [r,c] = obj.deviceTriplets(obj.map.genrouLin, ...
                @(i)obj.getRowIdx_fromGenrouIndices(i,0),nG);
            addSlice('genrouLin',r,c);

            % Nonlinear: harmonic-dense.
            Pnlin = { gi, gi      % d(Lsp*i*wr)/di
                      gi, gw      % d(Lsp*i*wr)/dwr
                      ge, ge      % d(invPark*v)/dv
                      ge, gt      % d(invPark*v)/d(deltaTheta)
                      gw, gi      % dTe/di
                      gia, gi(1:6) % d(invPark*i)/di
                      gia, gt };   % d(invPark*i)/d(deltaTheta)
            obj.map.genrouNlin = obj.localPatternHbDense(Pnlin,eqG,eqG);
            [r,c] = obj.deviceTriplets(obj.map.genrouNlin, ...
                @(i)obj.getRowIdx_fromGenrouIndices(i,0),nG);
            addSlice('genrouNlin',r,c);

            % Park-transform voltage reads the terminal bus voltage at the
            % SAME harmonic (a constant gain, so harmonic-diagonal).
            rp = []; cp = [];
            for k = 0:obj.K
                for g = 1:nG
                    base = obj.getRowIdx_fromGenrouIndices(g,k)-1;
                    vabc = obj.getRowIdxAbc_fromBusNum(obj.genrou_bus(g),k);
                    [rr,cc] = saltStampTriplets(base+ge,vabc);
                    rp = [rp;rr]; cp = [cp;cc];                  %#ok<AGROW>
                end
            end
            addSlice('genrouParkV',rp,cp);

            % Only the DC angle equation sees w0. There is no w0 term at
            % the higher harmonics by the definition of the fundamental.
            rs = []; cs = [];
            for g = 1:nG
                base = obj.getRowIdx_fromGenrouIndices(g,0)-1;
                rs = [rs; base+gt(:)];                           %#ok<AGROW>
                cs = [cs; repmat(idxW0,numel(gt),1)];            %#ok<AGROW>
            end
            addSlice('genrouW0',rs,cs);

            %% ---------------- EXDC1 -----------------------------------
            e1 = idxHbExcDc1_x1; e2 = idxHbExcDc1_x2; e3 = idxHbExcDc1_x3;
            e4 = idxHbExcDc1_x4; e5 = idxHbExcDc1_x5; ev = idxHbExcDc1_Vt;

            Plin = { e1, e1     % -Ke - j*wk*Te
                     e1, e4     % +x4
                     e2, e2     % 1 + j*wk*Tr
                     e2, ev     % -et
                     e3, e2     % 1 + j*wk*Tc
                     e3, e3     % 1 + j*wk*Tb
                     e3, e5     % 1 + j*wk*Tc
                     e4, e3     % -Ka
                     e4, e4     % 1 + j*wk*Ta
                     e5, e1     % -j*wk*Kf1
                     e5, e5 };  % 1 + j*wk*Tf1
            obj.map.excDc1Lin = obj.localPatternHbDiag(Plin,eqE);
            [r,c] = obj.deviceTriplets(obj.map.excDc1Lin, ...
                @(i)obj.getRowIdx_fromExcDc1Indices(i,0),nE);
            addSlice('excDc1Lin',r,c);

            % Terminal voltage is |vabc|, a time-domain nonlinearity, so
            % both its own column and the bus columns are harmonic-dense.
            obj.map.excDc1Nlin = obj.localPatternHbDense({ev,ev},eqE,eqE);
            [r,c] = obj.deviceTriplets(obj.map.excDc1Nlin, ...
                @(i)obj.getRowIdx_fromExcDc1Indices(i,0),nE);
            addSlice('excDc1Nlin',r,c);

            rV = []; cV = []; rF = []; cF = [];
            for j = 1:nE
                base  = obj.getRowIdx_fromExcDc1Indices(j,0)-1;
                gbase = obj.getRowIdx_fromGenrouIndices( ...
                    obj.excDc1_genrouIdx(j),0)-1;
                rowsVt = reshape((0:obj.K)*eqE + base + ev(:),[],1);
                vabc = reshape(obj.getRowIdxAbc_fromBusNum( ...
                    obj.excDc1_bus(j),0:obj.K),[],1);
                % getRowIdxAbc_fromBusNum over a row of k already returns
                % harmonic-outer / node-inner, which is the order the
                % kernel's dense block uses.
                [rr,cc] = saltStampTriplets(rowsVt,vabc);
                rV = [rV;rr]; cV = [cV;cc];                      %#ok<AGROW>
                % The exciter drives the machine's efd equation: three
                % dense blocks, in the kernel's order.
                rowsEfd = reshape((0:obj.K)*eqG + gbase + gE(:),[],1);
                colsX1  = reshape((0:obj.K)*eqE + base  + e1(:),[],1);
                colsEfd = rowsEfd;
                colsWr  = reshape((0:obj.K)*eqG + gbase + gw(:),[],1);
                for cols = {colsX1,colsEfd,colsWr}
                    [rr,cc] = saltStampTriplets(rowsEfd,cols{1});
                    rF = [rF;rr]; cF = [cF;cc];                  %#ok<AGROW>
                end
            end
            addSlice('excDc1Vabc',rV,cV);
            addSlice('excDc1Efd',rF,cF);

            %% ---------------- TGOV1D ----------------------------------
            t1 = idxHbTgov1d_x1; t2 = idxHbTgov1d_x2;
            tp = idxHbTgov1d_Pmech;

            Plin = { t1, t1     % 1 + j*wk*T3
                     t1, t2     % -1 - j*wk*T2
                     t2, t2     % -R - j*wk*R*T1
                     tp, t1     % +x1
                     tp, tp };  % -pmech
            obj.map.tgov1dLin = obj.localPatternHbDiag(Plin,eqV);
            [r,c] = obj.deviceTriplets(obj.map.tgov1dLin, ...
                @(i)obj.getRowIdx_fromTgov1dIndices(i,0),nV);
            addSlice('tgov1dLin',r,c);

            % Droop and damping read the machine speed at the same
            % harmonic; the torque handoff is a DFT sandwich and is dense.
            rW = []; cW = []; rT = []; cT = [];
            for j = 1:nV
                base  = obj.getRowIdx_fromTgov1dIndices(j,0)-1;
                gbase = obj.getRowIdx_fromGenrouIndices( ...
                    obj.tgov1d_genrouIdx(j),0)-1;
                for k = 0:obj.K
                    for locRows = {t2,tp}
                        [rr,cc] = saltStampTriplets( ...
                            base+k*eqV+locRows{1}, gbase+k*eqG+gw);
                        rW = [rW;rr]; cW = [cW;cc];              %#ok<AGROW>
                    end
                end
                rowsTm  = reshape((0:obj.K)*eqG + gbase + gT(:),[],1);
                colsPm  = reshape((0:obj.K)*eqV + base  + tp(:),[],1);
                colsWr  = reshape((0:obj.K)*eqG + gbase + gw(:),[],1);
                colsTm  = rowsTm;
                for cols = {colsPm,colsWr,colsTm}
                    [rr,cc] = saltStampTriplets(rowsTm,cols{1});
                    rT = [rT;rr]; cT = [cT;cc];                  %#ok<AGROW>
                end
            end
            addSlice('tgov1dWr',rW,cW);
            addSlice('tgov1dTm',rT,cT);

            %% ---------------- GFM IBR ---------------------------------
            % Column groups, straight from HBINITINDICES.
            bi  = idxHbIbrGFM_IphABC;  bvg = idxHbIbrGFM_Vg;
            bvo = idxHbIbrGFM_Vo;      bif = idxHbIbrGFM_If;
            bio = idxHbIbrGFM_Io;      bvi = idxHbIbrGFM_Vi;
            bga = idxHbIbrGFM_Gamma;   bir = idxHbIbrGFM_IfRef;
            bX  = idxHbIbrGFM_X;       bwr = idxHbIbrGFM_Wr;
            bvq = idxHbIbrGFM_VoRefDroop_q;
            % P and Q are two independent first-order filters, so the PqAvg
            % pair is addressed in halves.
            bpa = idxHbIbrGFM_PqAvg(1:2); bqa = idxHbIbrGFM_PqAvg(3:4);
            bth = idxHbIbrGFM_DeltaTheta;
            bva = idxHbIbrGFM_Vo_VR_avg; bvv = idxHbIbrGFM_Vo_VR;
            bvr = idxHbIbrGFM_VoRef;   bio0 = idxHbIbrGFM_Io0;
            % EQUATION rows are NOT the same permutation as the variables:
            % the io differential equation sits on rows 11:14 while the io
            % VARIABLE occupies columns 19:22. Numerically each equation
            % group lands on the slot of the variable at its own ordinal in
            % the layout, which is why these aliases are built from the same
            % index vectors - the row order is the one tabulated in
            % STAMPOFHBIBRGFMMATRICES_FD and must stay in step with it.
            qPkI = bi;  qPkV = bvg; qCp  = bvo; qCap = bif; qFlt = bio;
            qCC  = bvi; qCI  = bga; qVC  = bir; qVI  = bX;  qDW  = bwr;
            qDV  = bvq; qPav = bpa; qQav = bqa; qAng = bth; qVRa = bva;
            qVRv = bvv; qVRr = bvr;
            % The dq0 zero-sequence pair, appended past the 60 rows the
            % device had when it was a 2x3 interface. The EQUATION on these
            % rows is v0 = 0 (a delta connection on the network side), and
            % the VARIABLE in these columns is the zero-sequence current io0
            % the inverter draws to hold it. See the DELTA, NOT WYE note at
            % the top of CALCHBIBRGFMMATRICES_FD.
            qPkZ = bio0;

            Plin = { qPkI, bi       %  1 -eye(3), Park I abc side
                     qPkV, bvg      %  2 -eye(2), Park V dq side
                     qCp,  bio      %  3 coupling, d/d(io)
                     qCp,  bvo      %  4
                     qCp,  bvg      %  5
                     qCap, bvo      %  6 capacitor, restored vo'
                     qCap, bif      %  7 1/Cf plus the Rcap ESR
                     qCap, bio      %  8
                     qFlt, bif      %  9 filter
                     qFlt, bvi      % 10
                     qFlt, bvo      % 11
                     qCC,  bvi      % 12 current controller
                     qCC,  bga      % 13
                     qCC,  bvo      % 14
                     qCI,  bga      % 15 current integrator error
                     qCI,  bir      % 16
                     qCI,  bif      % 17
                     qVC,  bir      % 18 voltage controller
                     qVC,  bX       % 19
                     qVC,  bio      % 20
                     qVI,  bX       % 21 voltage integrator error
                     qVI,  bvr      % 22
                     qVI,  bvo      % 23
                     qDW,  bwr      % 24 P droop
                     qDW,  bpa      % 25
                     qDV,  bvq      % 26 Q droop
                     qDV,  bqa      % 27
                     qPav, bpa      % 28 measured P
                     qQav, bqa      % 29 measured Q
                     qAng, bth      % 30 frame angle
                     qAng, bwr      % 31
                     qVRa, bva      % 32 VR filter
                     qVRa, bio      % 33
                     qVRv, bvv      % 34 VR drop
                     qVRv, bva      % 35
                     qVRv, bio      % 36
                     qVRr, bvr      % 37 VR reference
                     qVRr, bvv      % 38
                     qVRr, bvq      % 39
                     qPkI, bio0 };  % 40 dq0: the ones(3,1) 0-column
            obj.map.ibrGFMLin = obj.localPatternHbDiag(Plin,eqI);
            [r,c] = obj.deviceTriplets(obj.map.ibrGFMLin, ...
                @(i)obj.getRowIdx_fromIbrGFMIndices(i,0),nI);
            addSlice('ibrGFMLin',r,c);

            Pnlin = { qPkI, bio     %  1 d(invT2*io)/d(io)
                      qPkI, bth     %  2 d(invT2*io)/d(deltaTheta)
                      qPkV, bth     %  3 d(T2*vabc)/d(deltaTheta)
                      qCp,  bio     %  4 d(j*w*io)/d(io)
                      qCp,  bwr     %  5 d(j*w*io)/d(wr)
                      qCap, bvo     %  6
                      qCap, bwr     %  7
                      qFlt, bif     %  8
                      qFlt, bwr     %  9
                      qCC,  bif     % 10 d(-j*wr*Lf*if)/d(if)
                      qCC,  bwr     % 11
                      qVC,  bvo     % 12 d(-j*wr*Cf*vo)/d(vo)
                      qVC,  bwr     % 13
                      qPav, bvg     % 14 d(p)/d(vg)
                      qPav, bio     % 15
                      qQav, bvg     % 16 d(q)/d(vg)
                      qQav, bio };  % 17
            obj.map.ibrGFMNlin = obj.localPatternHbDense(Pnlin,eqI,eqI);
            [r,c] = obj.deviceTriplets(obj.map.ibrGFMNlin, ...
                @(i)obj.getRowIdx_fromIbrGFMIndices(i,0),nI);
            addSlice('ibrGFMNlin',r,c);

            % Park V reads the terminal bus voltage through T2(theta), a DFT
            % sandwich, so unlike the GENROU's identity gain this coupling is
            % harmonic-DENSE. Columns are harmonic-outer / node-inner, which
            % is what getRowIdxAbc_fromBusNum over a ROW of k returns and
            % what the kernel's dense block assumes.
            rV = []; cV = []; rs = []; cs = [];
            for j = 1:nI
                base  = obj.getRowIdx_fromIbrGFMIndices(j,0)-1;
                rowsV = reshape((0:obj.K)*eqI + base + qPkV(:),[],1);
                vabc  = reshape(obj.getRowIdxAbc_fromBusNum( ...
                    obj.ibrGFM_bus(j),0:obj.K),[],1);
                [rr,cc] = saltStampTriplets(rowsV,vabc);
                rV = [rV;rr]; cV = [cV;cc];                      %#ok<AGROW>
                % Frequency column. Every restored derivative pulls on w0,
                % so unlike the GENROU - whose only w0 column is the DC
                % angle row - the whole block is stamped, at every harmonic.
                rowsAll = reshape(obj.getRowIdx_fromIbrGFMIndices( ...
                    j,0:obj.K)-1 + (1:eqI).',[],1);
                rs = [rs;rowsAll];                               %#ok<AGROW>
                cs = [cs;repmat(idxW0,numel(rowsAll),1)];        %#ok<AGROW>
            end
            addSlice('ibrGFMVabc',rV,cV);
            addSlice('ibrGFMW0',rs,cs);

            % Zero-sequence row of the dq0 FORWARD Park, v0 = 0, written
            % against the terminal bus nodes. Its coefficients are
            % (1/3)*busToIbrGFM_V on the three phases and carry no theta, so
            % unlike Park V above this coupling is harmonic-DIAGONAL and its
            % values are filled once, with the other constants, below.
            rZ = []; cZ = [];
            for k = 0:obj.K
                for g = 1:nI
                    vabc = obj.getRowIdxAbc_fromBusNum(obj.ibrGFM_bus(g),k);
                    rowsZ = obj.getRowIdx_fromIbrGFMIndices(g,k)-1 + qPkZ;
                    [rr,cc] = saltStampTriplets(rowsZ,vabc);
                    rZ = [rZ;rr]; cZ = [cZ;cc];                  %#ok<AGROW>
                end
            end
            addSlice('ibrGFMZeroSeq',rZ,cZ);

            %% ---------------- Composite load (IM + ZIP) ---------------
            % Here the equation rows DO coincide with the variable layout,
            % so one set of index vectors serves as both.
            lv  = idxHbLoadComp_Vs;    li  = idxHbLoadComp_IphABC;
            lws = idxHbLoadComp_WSlip; lth = idxHbLoadComp_DeltaThetaSys;
            lis = idxHbLoadComp_Is;    lir = idxHbLoadComp_Ir;
            lwr = idxHbLoadComp_Wr;    lpq = idxHbLoadComp_PQ;
            lvm = idxHbLoadComp_Vm;    liz = idxHbLoadComp_Izip;
            % The stator/rotor equation is one 4x4 complex block over
            % isr = [is; ir], so its rows and columns are the two index
            % vectors taken together; only the STATOR half of it sees the
            % terminal voltage, because vsr = [vs; 0; 0].
            lsr = [lis lir];

            Plin = { lv,  lv        %  1 -eye(4), Park V dq side
                     li,  li        %  2 -eye(6), Park I abc side
                     lws, lws       %  3 slip
                     lws, lwr       %  4
                     lth, lth       %  5 frame angle, pinned
                     lis, lv        %  6 -vsr, stator rows only
                     lsr, lsr       %  7 R + j*wk*L_PU + wMatSys*L_PU
                     lwr, lwr       %  8 restored swing inertia
                     lpq, lpq       %  9 -eye(4), ZIP power
                     lvm, lvm       % 10 -(1 + j*wk*Tau)
                     liz, liz };    % 11 -eye(4), ZIP current
            obj.map.loadCompLin = obj.localPatternHbDiag(Plin,eqL);
            [r,c] = obj.deviceTriplets(obj.map.loadCompLin, ...
                @(i)obj.getRowIdx_fromLoadCompIndices(i,0),nL);
            addSlice('loadCompLin',r,c);

            Pnlin = { lv,  lth      %  1 d(T2*vabc)/d(deltaThetaSys)
                      li,  lis      %  2 d(invT2*(is+izip))/d(is)
                      li,  liz      %  3 d(invT2*(is+izip))/d(izip)
                      li,  lth      %  4 d(invT2*(is+izip))/d(deltaTheta)
                      lsr, lsr      %  5 d(wMatSlip*L_PU*isr)/d(isr)
                      lsr, lws      %  6 d(wMatSlip*L_PU*isr)/d(wslip)
                      lwr, lsr      %  7 d(isr'*L_Te*isr)/d(isr)
                      lwr, lwr      %  8 d(-T0*(wr/w0)^m)/d(wr)
                      lpq, lvm      %  9 d(PQz*Vm^2 + PQi*Vm)/d(Vm)
                      lvm, lv       % 10 d(|vs|)/d(vs)
                      liz, lpq      % 11 d(iZipEval/Vm^2)/d(PQ)
                      liz, lv       % 12 d(iZipEval/Vm^2)/d(vs)
                      liz, lvm };   % 13 d(iZipEval/Vm^2)/d(Vm)
            obj.map.loadCompNlin = obj.localPatternHbDense(Pnlin,eqL,eqL);
            [r,c] = obj.deviceTriplets(obj.map.loadCompNlin, ...
                @(i)obj.getRowIdx_fromLoadCompIndices(i,0),nL);
            addSlice('loadCompNlin',r,c);

            rV = []; cV = []; rs = []; cs = [];
            for j = 1:nL
                base  = obj.getRowIdx_fromLoadCompIndices(j,0)-1;
                rowsV = reshape((0:obj.K)*eqL + base + lv(:),[],1);
                vabc  = reshape(obj.getRowIdxAbc_fromBusNum( ...
                    obj.loadComp_bus(j),0:obj.K),[],1);
                [rr,cc] = saltStampTriplets(rowsV,vabc);
                rV = [rV;rr]; cV = [cV;cc];                      %#ok<AGROW>
                rowsAll = reshape(obj.getRowIdx_fromLoadCompIndices( ...
                    j,0:obj.K)-1 + (1:eqL).',[],1);
                rs = [rs;rowsAll];                               %#ok<AGROW>
                cs = [cs;repmat(idxW0,numel(rowsAll),1)];        %#ok<AGROW>
            end
            addSlice('loadCompVabc',rV,cV);
            addSlice('loadCompW0',rs,cs);

            %% ---------------- Slack-angle (frequency) equation --------
            vSlack = obj.getRowIdxAbc_fromBusNum(obj.slack_busNum,1);
            addSlice('phase',[idxW0;idxW0],vSlack(1:2));

            %% ---------------- Imaginary-DC enforcement ----------------
            % A real signal has a REAL DC coefficient, so every "imaginary
            % part at k = 0" unknown is structurally zero and its equation
            % is identically zero - the Jacobian would be rank deficient.
            % GAMMA's second row is identically zero for the same reason,
            % so every nonlinear block already contributes nothing there.
            %
            % With triplets the legacy's in-place overwrite is not
            % available (triplets accumulate), so instead: every value
            % landing in one of these rows is zeroed at assembly, and this
            % slice supplies the 1 on the diagonal. Cleaner than the
            % legacy, which left a few benign cross terms behind in those
            % rows; the solution is the same either way, which is why
            % validateHb.m compares Jacobians on the complement of this
            % row set and checks these entries are zero separately.
            rowsDc = obj.imagDcRowIdx;
            obj.map.isImagDcRow = false(obj.getSysSizeAll,1);
            obj.map.isImagDcRow(rowsDc) = true;
            addSlice('imagDc',rowsDc,rowsDc);

            %% ---------------- Finalise --------------------------------
            obj.stampRow = vertcat(R{:});
            obj.stampCol = vertcat(C{:});
            obj.stampVal = zeros(n,1);
            obj.stampSlice = slice;

            % Which triplets get killed at assembly. Both pinned row sets
            % are cleared the same way, and each then has its own slice put
            % the equation back.
            obj.map.isPinnedRow = ...
                obj.map.isImagDcRow;
            obj.map.killVal = obj.map.isPinnedRow(obj.stampRow);
            obj.map.killVal(obj.stampSlice.imagDc) = false;

            % Entries that never change, filled once
            obj.stampVal(obj.stampSlice.injGenrou) = ...
                repmat(reshape(-eye(6),[],1),nG*K1,1);
            % The inverter produces, so its terminal current leaves the
            % device; the composite load consumes, so its current flows in.
            obj.stampVal(obj.stampSlice.injIbrGFM) = ...
                repmat(reshape(-eye(6),[],1),nI*K1,1);
            obj.stampVal(obj.stampSlice.injLoadComp) = ...
                repmat(reshape(eye(6),[],1),nL*K1,1);
            % v0 = (va+vb+vc)/3 = 0, real part against the three real nodes
            % and imaginary against the three imaginary ones. The slice was
            % built k-outer / inverter-inner, so the per-inverter block
            % varies fastest.
            if nI > 0
                [~,busIdxI] = ismembertol(obj.ibrGFM_bus,obj.bus);
                cZeroSeq = obj.busToIbrGFM_V(busIdxI,(1:nI).')/3;
                blkZ = reshape(repmat(eye(2),1,3),[],1).*cZeroSeq(:).';
                obj.stampVal(obj.stampSlice.ibrGFMZeroSeq) = ...
                    repmat(blkZ(:),K1,1);
            end
            obj.stampVal(obj.stampSlice.xfmrIsrc) = ...
                repmat(reshape(blkdiag([eye(2);-eye(2)],[eye(2);-eye(2)], ...
                [eye(2);-eye(2)]),[],1),numel(obj.map.xfmr.idxMna)*K1,1);
            obj.stampVal(obj.stampSlice.phase) = ...
                [obj.slack_imScalConst; -1];
            obj.stampVal(obj.stampSlice.imagDc) = 1;
        end

        function idx = imagDcRowIdx(obj)
            %IMAGDCROWIDX  Rows holding the imaginary part of a DC term.
            idx = (2:2:obj.getLinSysSize).';
            if obj.numGenrou > 0
                r = obj.getRowIdx_fromGenrouIndices((1:obj.numGenrou).',0)-1;
                idx = [idx; reshape(r.' + (2:2:obj.numGenrouEqns).',[],1)];
            end
            if obj.numExcDc1 > 0
                r = obj.getRowIdx_fromExcDc1Indices((1:obj.numExcDc1).',0)-1;
                idx = [idx; reshape(r.' + (2:2:obj.numExcDc1Eqns).',[],1)];
            end
            if obj.numTgov1d > 0
                r = obj.getRowIdx_fromTgov1dIndices((1:obj.numTgov1d).',0)-1;
                idx = [idx; reshape(r.' + (2:2:obj.numTgov1dEqns).',[],1)];
            end
            if obj.numIbrGFM > 0
                r = obj.getRowIdx_fromIbrGFMIndices((1:obj.numIbrGFM).',0)-1;
                idx = [idx; reshape(r.' + (2:2:obj.numIbrGFMEqns).',[],1)];
            end
            if obj.numLoadComp > 0
                r = obj.getRowIdx_fromLoadCompIndices((1:obj.numLoadComp).',0)-1;
                idx = [idx; reshape(r.' + (2:2:obj.numLoadCompEqns).',[],1)];
            end
            idx = sort(idx);
        end
    end

    methods (Hidden)
        %% ------------- Helpers used by BUILDSTAMPMAPS -----------------
        function idxFT = branchRowIdxFT(obj,busFrom,busTo,k)
            % Interleave from/to node rows phase by phase, so a branch
            % block reads [Fa_re Fa_im Ta_re Ta_im Fb_re ...].
            idxF = obj.getRowIdxAbc_fromBusNum(busFrom,k);
            idxT = obj.getRowIdxAbc_fromBusNum(busTo,k);
            idxFT = reshape([idxF(1:2:end),idxF(2:2:end), ...
                             idxT(1:2:end),idxT(2:2:end)].',[],1);
        end

        function m = localPatternHbDiag(obj,P,nEqns)
            % Harmonic-DIAGONAL local pattern: the block list is replicated
            % once per harmonic with no cross-harmonic entries. Flat order
            % is harmonic-outer, entry-inner, which is the order the
            % per-harmonic linear kernels are called in.
            row = []; col = [];
            for k = 0:obj.K
                for j = 1:size(P,1)
                    [rr,cc] = saltStampTriplets(k*nEqns+P{j,1}(:), ...
                                                k*nEqns+P{j,2}(:));
                    row = [row;rr]; col = [col;cc];              %#ok<AGROW>
                end
            end
            m.row = row; m.col = col; m.nVals = numel(row);
            m.S = sparse(row,1:numel(row),1,nEqns*(obj.K+1),numel(row));
        end

        function m = localPatternHbDense(obj,P,nEqnsRow,nEqnsCol)
            % Harmonic-DENSE local pattern: each block couples every
            % harmonic of its rows to every harmonic of its columns, which
            % is what GAMMA*diag(.)*invGAMMA produces. Rows and columns are
            % harmonic-outer / variable-inner so the kernel's dense block
            % can be flattened straight in.
            row = []; col = [];
            for j = 1:size(P,1)
                rH = reshape((0:obj.K)*nEqnsRow + P{j,1}(:),[],1);
                cH = reshape((0:obj.K)*nEqnsCol + P{j,2}(:),[],1);
                [rr,cc] = saltStampTriplets(rH,cH);
                row = [row;rr]; col = [col;cc];                  %#ok<AGROW>
            end
            m.row = row; m.col = col; m.nVals = numel(row);
            m.S = sparse(row,1:numel(row),1,nEqnsRow*(obj.K+1),numel(row));
        end

        function [r,c] = deviceTriplets(~,m,rowBaseFcn,nDev)
            if nDev==0, r = []; c = []; return, end
            base = rowBaseFcn((1:nDev).') - 1;
            r = reshape(m.row + base.',[],1);
            c = reshape(m.col + base.',[],1);
        end

        function [VF,VT] = branchVoltages(obj,busFrom,busTo,k)
            nb = numel(busFrom);
            iF = obj.getRowIdxAbc_fromBusNums(busFrom,k);
            iT = obj.getRowIdxAbc_fromBusNums(busTo,k);
            % getRowIdxAbc_fromBusNums is phase-major, so the gather runs
            % bus-fastest within each phase: reshape to (nb x 3) and
            % transpose, NOT straight to (3 x nb).
            VF = reshape(complex(obj.x_v(iF(1:2:end)), ...
                obj.x_v(iF(2:2:end))),nb,3).';
            VT = reshape(complex(obj.x_v(iT(1:2:end)), ...
                obj.x_v(iT(2:2:end))),nb,3).';
        end

        function Iout = accumI(~,Iout,rows,vals)
            Iout = Iout + accumarray(rows(:),vals(:),size(Iout));
        end

        function P = gatherPages(obj,rowBase,loc)
            X = obj.x_v(rowBase(:).' + loc(:));
            P = reshape(X,numel(loc),1,numel(rowBase));
        end

        function Jx = blockJx(obj,m,rowBase,Jvals)
            % J(x)*x for every device, straight from the triplets.
            xLoc = obj.x_v(m.col + rowBase(:).');
            Jx = m.S * (Jvals .* xLoc);
        end
    end

    methods (Static, Hidden)
        %% ------------- Real/imaginary block helpers -------------------
        % These mirror SALT's static helpers. They are duplicated rather
        % than called across classes so HB stands on its own.
        function B = pairBlock2x2(a11,a12,a21,a22)
            nb = numel(a11);
            B = zeros(6,6,nb);
            for ph = 0:2
                B(2*ph+1,2*ph+1,:) = a11;  B(2*ph+1,2*ph+2,:) = a12;
                B(2*ph+2,2*ph+1,:) = a21;  B(2*ph+2,2*ph+2,:) = a22;
            end
        end

        function X = interleaveFT(X6)
            X = X6([1 4 2 5 3 6],:);
        end

        function v = interleaveReIm(Xc)
            [p,nb] = size(Xc);
            v = zeros(2*p,nb);
            v(1:2:end,:) = real(Xc);
            v(2:2:end,:) = imag(Xc);
            v = v(:);
        end

        function E = reImDropLast(M)
            % p x q x n complex -> 2p x (2q-1) x n real. The last column is
            % dropped because the extra unknown (w0) is real, so only its
            % real-part column means anything.
            [p,q,nb] = size(M);
            E = zeros(2*p,2*q,nb);
            E(1:2:end,1:2:end,:) =  real(M);
            E(1:2:end,2:2:end,:) = -imag(M);
            E(2:2:end,1:2:end,:) =  imag(M);
            E(2:2:end,2:2:end,:) =  real(M);
            E(:,end,:) = [];
        end

        function [vals,Ivals] = networkBlock(dIdV,dIdw,Vstamp,Istamp,w0)
            nb = size(dIdV,3); p = size(dIdV,1);
            M = cat(2,dIdV,reshape(dIdw,p,1,nb));
            E = HB.reImDropLast(M);
            vals = E(:);
            Ic = pagemtimes(dIdV,reshape(Vstamp,p,1,nb)) ...
               + reshape(dIdw,p,1,nb)*w0 - reshape(Istamp,p,1,nb);
            Ivals = HB.interleaveReIm(reshape(Ic,p,nb));
        end

        function Yeq = calcYeq_RLSeries(R,L,wk), Yeq = 1./(R + 1j*wk*L); end
        function Yeq = calcYeq_C(C,wk),          Yeq = 1j*wk*C;          end
        function Yeq = calcYeq_RParC(R,C,wk),    Yeq = 1./R + 1j*wk*C;   end
        function Zeq = calcZeq_L(L,wk),          Zeq = 1j*wk*L;          end

        % d/dw0 of each element admittance. The factor k is the chain rule
        % through w_k = k*w0, and it is why a higher harmonic pulls harder
        % on the frequency unknown than the fundamental does.
        function dIdw0 = calc_dIdw0_RLSeries(R,L,V,wk,k)
            Z = R + 1j*wk*L;
            dIdw0 = -(1j*k*L)./(Z.^2).*V;
        end
        function dIdw0 = calc_dIdw0_C(C,V,k), dIdw0 = 1j*k*C.*V; end
    end

    methods
        %% ========================= Stamping ==========================
        function obj = stampLinear_Y_I(obj)
            %STAMPLINEAR_Y_I  The passive network, one harmonic at a time.
            %
            % The network is LTI, so harmonic k sees only harmonic k and
            % this part of the Jacobian is block diagonal. The exception is
            % the w0 column: because every tone is w_k = k*w0, each
            % harmonic's admittance derivative carries a factor k, so the
            % higher harmonics pull harder on the frequency unknown than
            % the fundamental does.
            obj.I = zeros(obj.getSysSizeAll,1);
            w0 = obj.x_v(obj.getRowIdx_w0);
            K1 = obj.K+1;

            nTx  = numel(obj.tx_busFrom);
            iA   = obj.map.xfmr.idxAdm;
            iM   = obj.map.xfmr.idxMna;
            nLd  = numel(obj.load_bus);

            valTx = zeros(12*13*nTx,K1);
            valXA = zeros(12*13*numel(iA),K1);
            valLd = zeros(6*7*nLd,K1);
            % The MNA transformer row carries 19 real columns, not 13: its
            % block is [dF/dV (6) | dF/dI (3) | dF/dw0 (1)] = 10 complex
            % columns, and reImDropLast turns that into 2*10-1.
            valXM = zeros(6*19*numel(iM),K1);

            for k = 0:obj.K
                wk = k*w0;

                % ---- Transmission lines (three-phase PI) -------------
                if nTx > 0
                    [VF,VT] = obj.branchVoltages(obj.tx_busFrom,obj.tx_busTo,k);
                    Yb  = (1./(obj.tx_R + 1j*wk*obj.tx_L)).';
                    Ysh = (1j*wk*obj.tx_C).';
                    Ibr = Yb.*(VF-VT);
                    Iabc = [Ibr + Ysh.*VF; -Ibr + Ysh.*VT];
                    dIdV = HB.pairBlock2x2(Yb+Ysh,-Yb,-Yb,Yb+Ysh);
                    dBr  = HB.calc_dIdw0_RLSeries(obj.tx_R.',obj.tx_L.', ...
                        VF-VT,wk,k);
                    dIdw = [dBr + HB.calc_dIdw0_C(obj.tx_C.',VF,k); ...
                           -dBr + HB.calc_dIdw0_C(obj.tx_C.',VT,k)];
                    [valTx(:,k+1),Iv] = HB.networkBlock(dIdV, ...
                        HB.interleaveFT(dIdw),HB.interleaveFT([VF;VT]), ...
                        HB.interleaveFT(Iabc),w0);
                    obj.I = obj.accumI(obj.I,obj.map.txIrows(:,k+1),Iv);
                end

                % ---- Transformers with resistance (admittance) --------
                if ~isempty(iA)
                    [VF,VT] = obj.branchVoltages(obj.xfmr_busFrom(iA), ...
                        obj.xfmr_busTo(iA),k);
                    Ym = (1./(obj.xfmr_R(iA) + 1j*wk*obj.xfmr_L(iA))).';
                    Iabc3 = Ym.*(VF-VT);
                    dIdV = HB.pairBlock2x2(Ym,-Ym,-Ym,Ym);
                    dIdw = HB.calc_dIdw0_RLSeries(obj.xfmr_R(iA).', ...
                        obj.xfmr_L(iA).',VF-VT,wk,k);
                    [valXA(:,k+1),Iv] = HB.networkBlock(dIdV, ...
                        HB.interleaveFT([dIdw;-dIdw]), ...
                        HB.interleaveFT([VF;VT]), ...
                        HB.interleaveFT([Iabc3;-Iabc3]),w0);
                    obj.I = obj.accumI(obj.I,obj.map.xfmrYIrows(:,k+1),Iv);
                end

                % ---- Lossless transformers (branch current unknown) ---
                if ~isempty(iM)
                    [VF,VT] = obj.branchVoltages(obj.xfmr_busFrom(iM), ...
                        obj.xfmr_busTo(iM),k);
                    nM = numel(iM);
                    L = obj.xfmr_L(iM).';
                    Z = 1j*wk*L;
                    Iabc = obj.mnaCurrents(iM,k);
                    % F = -(Vf - Vt) + Z*I. At k = 0 this collapses to
                    % Vf = Vt, which is right: an inductor is a DC short.
                    F = -(VF-VT) + Z.*Iabc;
                    dFdV = repmat(-blkdiag([1 -1],[1 -1],[1 -1]),1,1,nM);
                    dFdI = reshape(eye(3),3,3,1).*reshape(Z,1,1,nM);
                    % dZ/dw0 = j*k*L, so dF/dw0 = j*k*L*I
                    dFdw = reshape(1j*k*L.*Iabc,3,1,nM);
                    M = cat(2,dFdV,dFdI,dFdw);
                    E = HB.reImDropLast(M);
                    valXM(:,k+1) = E(:);
                    Ic = pagemtimes(dFdV,reshape( ...
                            HB.interleaveFT([VF;VT]),6,1,nM)) ...
                       + pagemtimes(dFdI,reshape(Iabc,3,1,nM)) ...
                       + dFdw*w0 - reshape(F,3,1,nM);
                    obj.I = obj.accumI(obj.I,obj.map.xfmrZIrows(:,k+1), ...
                        HB.interleaveReIm(reshape(Ic,3,nM)));
                end

                % ---- Constant impedances to ground -------------------
                if nLd > 0
                    idx = obj.getRowIdxAbc_fromBusNums(obj.load_bus,k);
                    % Phase-major gather: reshape to (nLd x 3) and
                    % transpose, not straight to (3 x nLd).
                    Vabc = reshape(complex(obj.x_v(idx(1:2:end)), ...
                        obj.x_v(idx(2:2:end))),nLd,3).';
                    [Ym,dYdw] = obj.constLoadAdmittance(wk,k);
                    Iabc = Ym.*Vabc;
                    dIdV = reshape(eye(3),3,3,1).*reshape(Ym,1,1,nLd);
                    dIdw = reshape(dYdw.*Vabc,3,1,nLd);
                    M = cat(2,dIdV,dIdw);
                    E = HB.reImDropLast(M);
                    valLd(:,k+1) = E(:);
                    Ic = pagemtimes(dIdV,reshape(Vabc,3,1,nLd)) ...
                       + dIdw*w0 - reshape(Iabc,3,1,nLd);
                    obj.I = obj.accumI(obj.I,obj.map.loadIrows(:,k+1), ...
                        HB.interleaveReIm(reshape(Ic,3,nLd)));
                end
            end

            if nTx > 0,      obj.stampVal(obj.stampSlice.txY)   = valTx(:); end
            if ~isempty(iA), obj.stampVal(obj.stampSlice.xfmrY) = valXA(:); end
            if ~isempty(iM), obj.stampVal(obj.stampSlice.xfmrZ) = valXM(:); end
            if nLd > 0,      obj.stampVal(obj.stampSlice.loadY) = valLd(:); end

            % ---- Harmonic current source --------------------------------
            % Right-hand side only: it is an independent source, so it has
            % no Jacobian entry. Without it the harmonics decouple from the
            % fundamental and the only solution is the harmonic-free one.
            if obj.includeHigherHarmonics
                for j = 1:numel(obj.higherHarmonics_bus)
                    rows = obj.getRowIdxAbc_fromBusNum( ...
                        obj.higherHarmonics_bus(j), ...
                        obj.higherHarmonics_harmonicNum(j));
                    obj.I(rows) = obj.I(rows) + obj.higherHarmonics_I(:,j);
                end
            end
        end

        function obj = stampPhase_Y_I(obj)
            %STAMPPHASE_Y_I  The equation that closes the autonomous problem.
            %
            % Adding w0 as an unknown adds one unknown, so it needs one
            % equation. A periodic solution of an autonomous system is
            % invariant under time shift, so without pinning the phase the
            % Jacobian is singular by exactly that one dimension. This
            % fixes the ANGLE of the slack bus's phase-A fundamental and
            % leaves its magnitude free:
            %
            %   0 = CONST*V_real + (-1)*V_imag
            %
            % The Jacobian entries are constant and were filled in
            % BUILDSTAMPMAPS; only the right-hand side is set here.
            obj.I(obj.getRowIdx_w0) = 0;
        end

        function obj = stampNonlinear_Y_I(obj)
            obj = obj.stampGenrou;
            obj = obj.stampExcDc1;
            obj = obj.stampTgov1d;
            obj = obj.stampIbrGFM;
            obj = obj.stampLoadComp;
        end

        function obj = assemble_Y_I(obj)
            %ASSEMBLE_Y_I  Sparse Jacobian from the (value,row,col) triplets.
            val = obj.stampVal;
            % Clear both structurally pinned row sets - the imaginary-DC
            % unknowns, and the DC zero-sequence reference when the case
            % needs one - and let each set's dedicated slice put its own
            % equation back. See the notes in BUILDSTAMPMAPS for why this
            % replaces the legacy's in-place overwrite.
            val(obj.map.killVal) = 0;
            obj.Y = saltAssembleY(obj.stampRow,obj.stampCol,val, ...
                obj.getSysSizeAll);
            % Both pinned equations have a zero right-hand side: an
            % imaginary DC coefficient is zero, and so is the slack bus's
            % DC zero-sequence voltage.
            obj.I(obj.map.isPinnedRow) = 0;
        end
    end

    methods (Hidden)
        function Iabc = mnaCurrents(obj,deviceIdx,k)
            nb = numel(deviceIdx);
            Iabc = zeros(3,nb);
            for j = 1:nb
                mna = obj.xfmr_mnaIdxMap( ...
                    obj.xfmr_mnaIdxMap(:,1)==deviceIdx(j),2);
                idxI = obj.getRowIdxAbc_MNA_fromMnaVarIdx(mna,k);
                Iabc(:,j) = complex(obj.x_v(idxI(1:2:end)), ...
                    obj.x_v(idxI(2:2:end)));
            end
        end

        function [Ym,dYdw] = constLoadAdmittance(obj,wk,k)
            % The nameplate reactance is defined at NOMINAL frequency, so
            % L and C are parameters derived once with w0nom; the harmonic
            % evaluation then uses w_k.
            nb = numel(obj.load_bus);
            Ym = zeros(1,nb); dYdw = zeros(1,nb);
            for j = 1:nb
                X = obj.load_X(j); Rr = obj.load_R(j);
                if X > 0
                    L = X/obj.w0nom;
                    Ym(j)   = HB.calcYeq_RLSeries(Rr,L,wk);
                    dYdw(j) = -(1j*k*L)/((Rr+1j*wk*L)^2);
                elseif X < 0
                    Rp = 1/real(1/obj.load_Z(j));
                    Cc = imag(1/obj.load_Z(j))/obj.w0nom;
                    Ym(j)   = HB.calcYeq_RParC(Rp,Cc,wk);
                    dYdw(j) = 1j*k*Cc;
                else
                    Ym(j)   = 1/Rr;
                    dYdw(j) = 0;
                end
            end
        end

        function obj = addBlockRhs(obj,rowBase,locRows,vals)
            rows = reshape(rowBase(:).' + locRows(:),[],1);
            obj.I(rows) = obj.I(rows) + vals(:);
        end

        function obj = writeBlockRhs(obj,rowBase,nEqns,vals)
            rows = reshape(rowBase(:).' + (1:nEqns).',[],1);
            obj.I(rows) = vals(:);
        end

        %% -------------------------- GENROU ----------------------------
        function obj = stampGenrou(obj)
            nG = obj.numGenrou;
            if nG==0, return, end
            global idxHbGenrou_Idq0 idxHbGenrou_Edq0 idxHbGenrou_Efd ...
                   idxHbGenrou_Wr idxHbGenrou_DeltaTheta ...
                   idxHbGenrou_IphABC idxHbGenrou_Tm                %#ok<GVMIS>

            K1 = obj.K+1;
            eqG = obj.numGenrouEqns;
            r = obj.getRowIdx_fromGenrouIndices((1:nG).',0)-1;
            w0 = obj.x_v(obj.getRowIdx_w0);
            [~,busIdx] = ismembertol(obj.genrou_bus,obj.bus);
            toGenV = reshape(obj.busToGen_V(busIdx,(1:nG).')* ...
                ones(nG,1),1,1,nG);
            toGenI = reshape(obj.busToGen_I(busIdx,(1:nG).'),1,1,nG);

            % ---- Linear part: once per harmonic --------------------
            nValsLin = obj.map.genrouLin.nVals/K1;
            JLin = zeros(nValsLin,K1,nG);
            Flin = zeros(eqG,K1,nG);
            dParkV = zeros(36,nG,K1);
            dThW0  = zeros(2,nG);
            for k = 0:obj.K
                rk = obj.getRowIdx_fromGenrouIndices((1:nG).',k)-1;
                [JLin(:,k+1,:),Fk,dPk,dWk] = stampOfHbGenrouMatrices_FD( ...
                    obj.gatherPages(rk,idxHbGenrou_Idq0), ...
                    obj.gatherPages(rk,idxHbGenrou_Edq0), ...
                    obj.gatherPages(rk,idxHbGenrou_Efd), ...
                    obj.gatherPages(rk,idxHbGenrou_Wr), ...
                    obj.gatherPages(rk,idxHbGenrou_DeltaTheta), ...
                    obj.gatherPages(rk,idxHbGenrou_IphABC), ...
                    obj.gatherPages(rk,idxHbGenrou_Tm), ...
                    k,w0, ...
                    obj.genrou_RMatrix,obj.genrou_LMatrix, ...
                    reshape(obj.genrou_H,1,1,nG), ...
                    reshape(obj.genrou_KD,1,1,nG), ...
                    obj.w0nom,obj.N,toGenV,toGenI);
                Flin(:,k+1,:) = reshape(Fk,eqG,1,nG);
                dParkV(:,:,k+1) = dPk;
                if k==0, dThW0 = dWk; end
            end
            obj.stampVal(obj.stampSlice.genrouLin)   = JLin(:);
            obj.stampVal(obj.stampSlice.genrouParkV) = dParkV(:);
            obj.stampVal(obj.stampSlice.genrouW0)    = dThW0(:);

            % ---- Nonlinear part: dense across harmonics -------------
            % devicePagesTd handles the (re,im) interleaving hbIdft needs
            % and returns [nVar x N x nDev] pages.
            gidx = @(loc) @(k) obj.getRowIdx_fromGenrouIndices( ...
                (1:nG).',k)-1 + loc;
            iTd  = obj.devicePagesTd(gidx(idxHbGenrou_Idq0),7,nG);
            vTd  = obj.devicePagesTd( ...
                gidx([idxHbGenrou_Edq0 idxHbGenrou_Efd]),4,nG);
            wrTd = obj.devicePagesTd(gidx(idxHbGenrou_Wr),1,nG);
            thTd = obj.devicePagesTd(gidx(idxHbGenrou_DeltaTheta),1,nG);

            [JNlin,Fnlin] = calcHbGenrouMatrices_FD( ...
                iTd,vTd,wrTd,thTd,w0, ...
                obj.genrou_LspMatrix,obj.w0nom, ...
                obj.K,obj.N,obj.GAMMA,obj.invGAMMA,obj.tSamp);
            obj.stampVal(obj.stampSlice.genrouNlin) = JNlin(:);

            % ---- Right-hand side ------------------------------------
            % For a LINEAR equation the Newton right-hand side
            % J*x - F(x) collapses to the equation's CONSTANT term, which
            % is exactly what stampOfHbGenrouMatrices_FD returns in Flin
            % (the -(KD/w)*w term of the swing equation, at DC only). So
            % the linear J*x must NOT be added here - and neither must the
            % terminal-voltage column or the w0 column, which are linear
            % too and cancel in the same way. Only the nonlinear terms
            % need the full J(x)*x - F(x). This is why the legacy adds
            % nothing to Igen for its Park-voltage and speed stamps.
            Jx_nlin = obj.blockJx(obj.map.genrouNlin,r,JNlin);
            obj = obj.writeBlockRhs(r,eqG*K1, ...
                Jx_nlin - Fnlin + reshape(Flin,eqG*K1,nG));
            % Rows 21:22 (efd) and 33:34 (Tm) are deliberately left at
            % zero: the exciter and governor own them and overwrite both
            % the Jacobian and the right-hand side there.
        end
    end

    methods (Hidden)
        %% -------------------------- EXDC1 -----------------------------
        function obj = stampExcDc1(obj)
            nE = obj.numExcDc1;
            if nE==0, return, end
            global idxHbExcDc1_x1 idxHbExcDc1_x2 idxHbExcDc1_x3 ...
                   idxHbExcDc1_x4 idxHbExcDc1_x5 idxHbExcDc1_Vt     %#ok<GVMIS>
            global idxHbGenrou_Efd idxHbGenrou_Wr                   %#ok<GVMIS>

            K1  = obj.K+1;
            eqE = obj.numExcDc1Eqns;
            w0  = obj.x_v(obj.getRowIdx_w0);
            r   = obj.getRowIdx_fromExcDc1Indices((1:nE).',0)-1;
            rg  = obj.getRowIdx_fromGenrouIndices(obj.excDc1_genrouIdx,0)-1;
            [~,busIdx] = ismembertol(obj.excDc1_bus,obj.bus);
            toGenV = reshape(obj.busToGen_V(busIdx,(1:nE).')* ...
                ones(nE,1),1,1,nE);

            % ---- Linear part: once per harmonic --------------------
            nValsLin = obj.map.excDc1Lin.nVals/K1;
            JLin = zeros(nValsLin,K1,nE);
            Flin = zeros(eqE,K1,nE);
            p3 = @(v) reshape(v,1,1,nE);
            for k = 0:obj.K
                rk = obj.getRowIdx_fromExcDc1Indices((1:nE).',k)-1;
                [JLin(:,k+1,:),Fk] = stampOfHbExcDc1Matrices_FD( ...
                    obj.gatherPages(rk,idxHbExcDc1_x1), ...
                    obj.gatherPages(rk,idxHbExcDc1_x2), ...
                    obj.gatherPages(rk,idxHbExcDc1_x3), ...
                    obj.gatherPages(rk,idxHbExcDc1_x4), ...
                    obj.gatherPages(rk,idxHbExcDc1_x5), ...
                    obj.gatherPages(rk,idxHbExcDc1_Vt), ...
                    k,w0, ...
                    p3(obj.excDc1_Ke),p3(obj.excDc1_Te),p3(obj.excDc1_Tr), ...
                    p3(obj.excDc1_Tc),p3(obj.excDc1_Tb),p3(obj.excDc1_Ta), ...
                    p3(obj.excDc1_Ka),p3(obj.excDc1_Tf1),p3(obj.excDc1_Kf1), ...
                    p3(obj.excDc1_vref0));
                Flin(:,k+1,:) = reshape(Fk,eqE,1,nE);
            end
            obj.stampVal(obj.stampSlice.excDc1Lin) = JLin(:);

            % ---- Nonlinear part -----------------------------------
            x1Td = obj.devicePagesTd(@(k)obj.getRowIdx_fromExcDc1Indices( ...
                (1:nE).',k)-1 + idxHbExcDc1_x1, 1,nE);
            VtTd = obj.devicePagesTd(@(k)obj.getRowIdx_fromExcDc1Indices( ...
                (1:nE).',k)-1 + idxHbExcDc1_Vt, 1,nE);
            vabcTd = obj.busPagesTd(obj.excDc1_bus,nE);
            efdTd = obj.devicePagesTd(@(k)obj.getRowIdx_fromGenrouIndices( ...
                obj.excDc1_genrouIdx,k)-1 + idxHbGenrou_Efd, 1,nE);
            wrTd = obj.devicePagesTd(@(k)obj.getRowIdx_fromGenrouIndices( ...
                obj.excDc1_genrouIdx,k)-1 + idxHbGenrou_Wr, 1,nE);

            [JNlin,Fnlin,JEfd,FEfd] = calcHbExcDc1Matrices_FD( ...
                x1Td,VtTd,vabcTd,efdTd,wrTd, ...
                obj.w0nom,reshape(obj.excDc1_genToExc,1,1,nE),toGenV, ...
                obj.K,obj.GAMMA,obj.invGAMMA);

            % The kernel returns the vt row's own column and the bus
            % columns as one array; they live in different slices because
            % the bus columns are not local to the device block.
            nVt   = (2*K1)^2;
            obj.stampVal(obj.stampSlice.excDc1Nlin) = JNlin(1:nVt,:);
            obj.stampVal(obj.stampSlice.excDc1Vabc) = JNlin(nVt+1:end,:);
            obj.stampVal(obj.stampSlice.excDc1Efd)  = JEfd(:);

            % ---- Right-hand side ----------------------------------
            % Rows 1:10 are linear, so their right-hand side is the
            % constant term alone (vref at DC). The vt row is nonlinear
            % and needs J(x)*x - F(x) over both its own column and the
            % bus columns.
            obj = obj.writeBlockRhs(r,eqE*K1,reshape(Flin,eqE*K1,nE));

            Vt   = obj.harmonicStack(r,idxHbExcDc1_Vt,eqE,nE);
            Vabc = obj.harmonicStackBus(obj.excDc1_bus,nE);
            Jvt   = reshape(JNlin(1:nVt,:),2*K1,2*K1,nE);
            Jvabc = reshape(JNlin(nVt+1:end,:),2*K1,6*K1,nE);
            rhsVt = pagemtimes(Jvt,Vt) + pagemtimes(Jvabc,Vabc) ...
                  - reshape(Fnlin(obj.harmonicLocalRows(idxHbExcDc1_Vt,eqE),:), ...
                            2*K1,1,nE);
            obj = obj.writeHarmonicRows(r,idxHbExcDc1_Vt,eqE, ...
                reshape(rhsVt,2*K1,nE));

            % The machine's field-voltage row belongs to the exciter, so
            % OVERWRITE it: the GENROU kernel left it empty on purpose.
            X1  = obj.harmonicStack(r, idxHbExcDc1_x1,eqE,nE);
            Efd = obj.harmonicStack(rg,idxHbGenrou_Efd,obj.numGenrouEqns,nE);
            Wr  = obj.harmonicStack(rg,idxHbGenrou_Wr, obj.numGenrouEqns,nE);
            Jb  = reshape(JEfd,2*K1,2*K1,3,nE);
            rhsEfd = pagemtimes(reshape(Jb(:,:,1,:),2*K1,2*K1,nE),X1) ...
                   + pagemtimes(reshape(Jb(:,:,2,:),2*K1,2*K1,nE),Efd) ...
                   + pagemtimes(reshape(Jb(:,:,3,:),2*K1,2*K1,nE),Wr) ...
                   - reshape(FEfd,2*K1,1,nE);
            obj = obj.writeHarmonicRows(rg,idxHbGenrou_Efd, ...
                obj.numGenrouEqns,reshape(rhsEfd,2*K1,nE));
        end

        %% -------------------------- TGOV1D ----------------------------
        function obj = stampTgov1d(obj)
            nV = obj.numTgov1d;
            if nV==0, return, end
            global idxHbTgov1d_x1 idxHbTgov1d_x2 idxHbTgov1d_Pmech  %#ok<GVMIS>
            global idxHbGenrou_Wr idxHbGenrou_Tm                    %#ok<GVMIS>

            K1  = obj.K+1;
            eqV = obj.numTgov1dEqns;
            eqG = obj.numGenrouEqns;
            w0  = obj.x_v(obj.getRowIdx_w0);
            r   = obj.getRowIdx_fromTgov1dIndices((1:nV).',0)-1;
            rg  = obj.getRowIdx_fromGenrouIndices(obj.tgov1d_genrouIdx,0)-1;
            p3  = @(v) reshape(v,1,1,nV);

            % ---- Linear part: once per harmonic --------------------
            nValsLin = obj.map.tgov1dLin.nVals/K1;
            JLin = zeros(nValsLin,K1,nV);
            Flin = zeros(eqV,K1,nV);
            dWr  = zeros(4,2,K1,nV);     % [x2 row; pmech row] couplings
            for k = 0:obj.K
                rk = obj.getRowIdx_fromTgov1dIndices((1:nV).',k)-1;
                [JLin(:,k+1,:),Fk,dFx2,dFpm] = stampOfHbTgov1dMatrices_FD( ...
                    obj.gatherPages(rk,idxHbTgov1d_x1), ...
                    obj.gatherPages(rk,idxHbTgov1d_x2), ...
                    obj.gatherPages(rk,idxHbTgov1d_Pmech), ...
                    k,w0, ...
                    p3(obj.tgov1d_T1),p3(obj.tgov1d_T2),p3(obj.tgov1d_T3), ...
                    p3(obj.tgov1d_R),p3(obj.tgov1d_Dt), ...
                    obj.w0nom,obj.N,p3(obj.tgov1d_Pref0));
                Flin(:,k+1,:) = reshape(Fk,eqV,1,nV);
                dWr(:,1,k+1,:) = reshape(dFx2,4,1,1,nV);
                dWr(:,2,k+1,:) = reshape(dFpm,4,1,1,nV);
            end
            obj.stampVal(obj.stampSlice.tgov1dLin) = JLin(:);
            % tgov1dWr was built device-outer, then harmonic, then
            % {x2 row, pmech row}, then the 4 entries of each 2x2 block -
            % which is exactly dWr's own dimension order.
            obj.stampVal(obj.stampSlice.tgov1dWr) = dWr(:);

            % ---- Nonlinear part -----------------------------------
            pmechTd = obj.devicePagesTd(@(k)obj.getRowIdx_fromTgov1dIndices( ...
                (1:nV).',k)-1 + idxHbTgov1d_Pmech, 1,nV);
            TmTd = obj.devicePagesTd(@(k)obj.getRowIdx_fromGenrouIndices( ...
                obj.tgov1d_genrouIdx,k)-1 + idxHbGenrou_Tm, 1,nV);
            wrTd = obj.devicePagesTd(@(k)obj.getRowIdx_fromGenrouIndices( ...
                obj.tgov1d_genrouIdx,k)-1 + idxHbGenrou_Wr, 1,nV);

            [JTm,FTm] = calcHbTgov1dMatrices_FD(pmechTd,TmTd,wrTd, ...
                obj.w0nom,obj.genP,obj.K,obj.GAMMA,obj.invGAMMA);
            obj.stampVal(obj.stampSlice.tgov1dTm) = JTm(:);

            % ---- Right-hand side ----------------------------------
            % All three governor rows are linear (the machine-speed column
            % cancels), so their right-hand side is the constant alone.
            obj = obj.writeBlockRhs(r,eqV*K1,reshape(Flin,eqV*K1,nV));

            % The machine's mechanical-torque row belongs to the governor.
            Pm = obj.harmonicStack(r, idxHbTgov1d_Pmech,eqV,nV);
            Wr = obj.harmonicStack(rg,idxHbGenrou_Wr,   eqG,nV);
            Tm = obj.harmonicStack(rg,idxHbGenrou_Tm,   eqG,nV);
            Jb = reshape(JTm,2*K1,2*K1,3,nV);
            rhsTm = pagemtimes(reshape(Jb(:,:,1,:),2*K1,2*K1,nV),Pm) ...
                  + pagemtimes(reshape(Jb(:,:,2,:),2*K1,2*K1,nV),Wr) ...
                  + pagemtimes(reshape(Jb(:,:,3,:),2*K1,2*K1,nV),Tm) ...
                  - reshape(FTm,2*K1,1,nV);
            obj = obj.writeHarmonicRows(rg,idxHbGenrou_Tm,eqG, ...
                reshape(rhsTm,2*K1,nV));
        end

        %% ------------------------- GFM IBR ----------------------------
        function obj = stampIbrGFM(obj)
            nI = obj.numIbrGFM;
            if nI==0, return, end
            global idxHbIbrGFM_IphABC idxHbIbrGFM_Vg idxHbIbrGFM_Vo ...
                   idxHbIbrGFM_If idxHbIbrGFM_Io idxHbIbrGFM_Vi ...
                   idxHbIbrGFM_Gamma idxHbIbrGFM_IfRef idxHbIbrGFM_X ...
                   idxHbIbrGFM_Wr idxHbIbrGFM_VoRefDroop_q ...
                   idxHbIbrGFM_PqAvg idxHbIbrGFM_DeltaTheta ...
                   idxHbIbrGFM_Vo_VR_avg idxHbIbrGFM_Vo_VR ...
                   idxHbIbrGFM_VoRef idxHbIbrGFM_Io0                %#ok<GVMIS>

            K1  = obj.K+1;
            eqI = obj.numIbrGFMEqns;
            w0  = obj.x_v(obj.getRowIdx_w0);
            r   = obj.getRowIdx_fromIbrGFMIndices((1:nI).',0)-1;
            [~,busIdx] = ismembertol(obj.ibrGFM_bus,obj.bus);
            p3 = @(v) reshape(v,1,1,nI);
            toIbrV = p3(obj.busToIbrGFM_V(busIdx,(1:nI).'));
            toBusI = p3(obj.ibrGFMToBus_I(busIdx,(1:nI).'));

            % ---- Linear part: once per harmonic --------------------
            nValsLin = obj.map.ibrGFMLin.nVals/K1;
            JLin = zeros(nValsLin,K1,nI);
            Flin = zeros(eqI,K1,nI);
            dW0  = zeros(eqI,K1,nI);   % the half that reaches the RHS
            dW0b = zeros(eqI,K1,nI);   % the bare-additive half, J only
            for k = 0:obj.K
                rk = obj.getRowIdx_fromIbrGFMIndices((1:nI).',k)-1;
                [JLin(:,k+1,:),Fk,dWk,dWkBare] = ...
                    stampOfHbIbrGFMMatrices_FD( ...
                    obj.gatherPages(rk,idxHbIbrGFM_IphABC), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Vg), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Vo), ...
                    obj.gatherPages(rk,idxHbIbrGFM_If), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Io), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Vi), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Gamma), ...
                    obj.gatherPages(rk,idxHbIbrGFM_IfRef), ...
                    obj.gatherPages(rk,idxHbIbrGFM_X), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Wr), ...
                    obj.gatherPages(rk,idxHbIbrGFM_VoRefDroop_q), ...
                    obj.gatherPages(rk,idxHbIbrGFM_PqAvg), ...
                    obj.gatherPages(rk,idxHbIbrGFM_DeltaTheta), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Vo_VR_avg), ...
                    obj.gatherPages(rk,idxHbIbrGFM_Vo_VR), ...
                    obj.gatherPages(rk,idxHbIbrGFM_VoRef), ...
                    k,w0, ...
                    p3(obj.ibrGFM_Lf),p3(obj.ibrGFM_Rf), ...
                    p3(obj.ibrGFM_Cf),p3(obj.ibrGFM_Rcap), ...
                    p3(obj.ibrGFM_Lc),p3(obj.ibrGFM_Rc), ...
                    p3(obj.ibrGFM_kC_i),p3(obj.ibrGFM_kC_p), ...
                    p3(obj.ibrGFM_GC), ...
                    p3(obj.ibrGFM_kV_i),p3(obj.ibrGFM_kV_p), ...
                    p3(obj.ibrGFM_GV), ...
                    p3(obj.ibrGFM_Mp),p3(obj.ibrGFM_Mq), ...
                    p3(obj.ibrGFM_pset),p3(obj.ibrGFM_qset), ...
                    p3(obj.ibrGFM_voDroop_q_set), ...
                    p3(obj.ibrGFM_wmeas),p3(obj.ibrGFM_w_VR), ...
                    p3(obj.ibrGFM_Rv), ...
                    toBusI, ...
                    obj.w0nom,obj.N);
                Flin(:,k+1,:) = reshape(Fk,eqI,1,nI);
                dW0(:,k+1,:)  = reshape(dWk,eqI,1,nI);
                dW0b(:,k+1,:) = reshape(dWkBare,eqI,1,nI);
            end
            obj.stampVal(obj.stampSlice.ibrGFMLin) = JLin(:);
            % The Jacobian's frequency column is the SUM of the two halves;
            % only dW0 goes on to the right-hand side below.
            obj.stampVal(obj.stampSlice.ibrGFMW0) = dW0(:) + dW0b(:);

            % ---- Nonlinear part: dense across harmonics -------------
            bidx = @(loc) @(k) obj.getRowIdx_fromIbrGFMIndices( ...
                (1:nI).',k)-1 + loc;
            vgTd = obj.devicePagesTd(bidx(idxHbIbrGFM_Vg),2,nI);
            voTd = obj.devicePagesTd(bidx(idxHbIbrGFM_Vo),2,nI);
            ifTd = obj.devicePagesTd(bidx(idxHbIbrGFM_If),2,nI);
            ioTd = obj.devicePagesTd(bidx(idxHbIbrGFM_Io),2,nI);
            wrTd = obj.devicePagesTd(bidx(idxHbIbrGFM_Wr),1,nI);
            thTd = obj.devicePagesTd(bidx(idxHbIbrGFM_DeltaTheta),1,nI);
            vabcTd = obj.busPagesTd(obj.ibrGFM_bus,nI);

            [JNlin,Fnlin,JVabc] = calcHbIbrGFMMatrices_FD( ...
                vgTd,voTd,ifTd,ioTd,wrTd,thTd,vabcTd,w0, ...
                p3(obj.ibrGFM_Lf),p3(obj.ibrGFM_Cf),obj.w0nom, ...
                toBusI,toIbrV, ...
                obj.K,obj.N,obj.GAMMA,obj.invGAMMA,obj.tSamp);
            obj.stampVal(obj.stampSlice.ibrGFMNlin) = JNlin(:);
            obj.stampVal(obj.stampSlice.ibrGFMVabc) = JVabc;

            % ---- Right-hand side ------------------------------------
            % Same accounting as the GENROU: a purely linear equation
            % contributes only its constant term, because J*x and F cancel
            % there. Two things break that cancellation here, so their J*x
            % IS added:
            %   the nonlinear block, as usual;
            %   the frequency column's state-dependent half dW0, which is
            %       the w0 side of terms BILINEAR in w0 and a state
            %       (j*k*w0*x). For those J*x is twice the term, not once,
            %       so dW0*w0 has to appear. dW0b - the bare additive -w0
            %       of the angle equation - is linear and must NOT.
            % The Park-V coupling into the bus voltage is nonlinear too, so
            % its J*x is added on the Park-V rows, exactly as the exciter
            % does for its terminal-voltage row.
            Jx = obj.blockJx(obj.map.ibrGFMNlin,r,JNlin);
            obj = obj.writeBlockRhs(r,eqI*K1, ...
                Jx - Fnlin + reshape(Flin,eqI*K1,nI) ...
                   + reshape(dW0,eqI*K1,nI)*w0);

            Vabc = obj.harmonicStackBus(obj.ibrGFM_bus,nI);
            JVb  = reshape(JVabc,4*K1,6*K1,nI);
            obj = obj.addBlockRhs(r, ...
                obj.harmonicLocalRows(idxHbIbrGFM_Vg,eqI), ...
                reshape(pagemtimes(JVb,Vabc),4*K1,nI));
        end

        %% ------------------ Composite load (IM + ZIP) -----------------
        function obj = stampLoadComp(obj)
            nL = obj.numLoadComp;
            if nL==0, return, end
            global idxHbLoadComp_Vs idxHbLoadComp_IphABC ...
                   idxHbLoadComp_WSlip idxHbLoadComp_DeltaThetaSys ...
                   idxHbLoadComp_Is idxHbLoadComp_Ir idxHbLoadComp_Wr ...
                   idxHbLoadComp_PQ idxHbLoadComp_Vm ...
                   idxHbLoadComp_Izip                               %#ok<GVMIS>

            K1  = obj.K+1;
            eqL = obj.numLoadCompEqns;
            w0  = obj.x_v(obj.getRowIdx_w0);
            r   = obj.getRowIdx_fromLoadCompIndices((1:nL).',0)-1;
            [~,busIdx] = ismembertol(obj.loadComp_bus,obj.bus);
            p3 = @(v) reshape(v,1,1,nL);
            % The voltage base change is identically 1 (see
            % BUSTOLOADCOMP_V), so it is expanded to one page per load the
            % same way the GENROU expands busToGen_V.
            toLoadV = reshape(obj.busToLoadComp_V(busIdx,(1:nL).')* ...
                ones(nL,1),1,1,nL);
            toBusI  = p3(obj.loadCompToBus_I(busIdx,(1:nL).'));
            % L_PU carries the TD model's normalisation, L/w0nom, so every
            % restored derivative is j*k*w0 times that same constant.
            L_PU = obj.loadCompIM_L/obj.w0nom;
            PQz = reshape(obj.loadCompZip_PQz,2,1,nL);
            PQi = reshape(obj.loadCompZip_PQi,2,1,nL);
            PQp = reshape(obj.loadCompZip_PQp,2,1,nL);

            % ---- Linear part: once per harmonic --------------------
            nValsLin = obj.map.loadCompLin.nVals/K1;
            JLin = zeros(nValsLin,K1,nL);
            Flin = zeros(eqL,K1,nL);
            dW0  = zeros(eqL,K1,nL);
            dW0b = zeros(eqL,K1,nL);
            for k = 0:obj.K
                rk = obj.getRowIdx_fromLoadCompIndices((1:nL).',k)-1;
                [JLin(:,k+1,:),Fk,~,dWk,dWkBare] = ...
                    stampOfHbLoadCompMatrices_FD( ...
                    obj.gatherPages(rk,idxHbLoadComp_Vs), ...
                    obj.gatherPages(rk,idxHbLoadComp_IphABC), ...
                    obj.gatherPages(rk,idxHbLoadComp_WSlip), ...
                    obj.gatherPages(rk,idxHbLoadComp_DeltaThetaSys), ...
                    obj.gatherPages(rk,idxHbLoadComp_Is), ...
                    obj.gatherPages(rk,idxHbLoadComp_Ir), ...
                    obj.gatherPages(rk,idxHbLoadComp_Wr), ...
                    obj.gatherPages(rk,idxHbLoadComp_PQ), ...
                    obj.gatherPages(rk,idxHbLoadComp_Vm), ...
                    obj.gatherPages(rk,idxHbLoadComp_Izip), ...
                    k,w0, ...
                    obj.loadCompIM_R,L_PU,obj.loadCompIM_LTe, ...
                    p3(obj.loadCompIM_Tm0),p3(obj.loadCompIM_m), ...
                    p3(obj.loadCompZip_Tau), PQz,PQi,PQp, ...
                    obj.w0nom,obj.N,p3(obj.loadCompIM_H), ...
                    toLoadV,toBusI);
                Flin(:,k+1,:) = reshape(Fk,eqL,1,nL);
                dW0(:,k+1,:)  = reshape(dWk,eqL,1,nL);
                dW0b(:,k+1,:) = reshape(dWkBare,eqL,1,nL);
            end
            obj.stampVal(obj.stampSlice.loadCompLin) = JLin(:);

            % ---- Nonlinear part: dense across harmonics -------------
            lidx = @(loc) @(k) obj.getRowIdx_fromLoadCompIndices( ...
                (1:nL).',k)-1 + loc;
            vsTd  = obj.devicePagesTd(lidx(idxHbLoadComp_Vs),2,nL);
            iTd   = obj.devicePagesTd(lidx(idxHbLoadComp_IphABC),3,nL);
            wslTd = obj.devicePagesTd(lidx(idxHbLoadComp_WSlip),1,nL);
            thTd  = obj.devicePagesTd( ...
                lidx(idxHbLoadComp_DeltaThetaSys),1,nL);
            isTd  = obj.devicePagesTd(lidx(idxHbLoadComp_Is),2,nL);
            irTd  = obj.devicePagesTd(lidx(idxHbLoadComp_Ir),2,nL);
            wrTd  = obj.devicePagesTd(lidx(idxHbLoadComp_Wr),1,nL);
            PQTd  = obj.devicePagesTd(lidx(idxHbLoadComp_PQ),2,nL);
            VmTd  = obj.devicePagesTd(lidx(idxHbLoadComp_Vm),1,nL);
            izTd  = obj.devicePagesTd(lidx(idxHbLoadComp_Izip),2,nL);
            vabcTd = obj.busPagesTd(obj.loadComp_bus,nL);

            [JNlin,Fnlin,JVabc,~,dW0nlin] = calcHbLoadCompMatrices_FD( ...
                vsTd,iTd,wslTd,thTd,isTd,irTd,wrTd,PQTd,VmTd,izTd, ...
                vabcTd,w0, ...
                L_PU,obj.loadCompIM_LTe, ...
                p3(obj.loadCompIM_Tm0),p3(obj.loadCompIM_m), ...
                PQz,PQi,PQp, toLoadV,toBusI, ...
                obj.K,obj.N,obj.GAMMA,obj.invGAMMA,obj.tSamp);
            obj.stampVal(obj.stampSlice.loadCompNlin) = JNlin(:);
            obj.stampVal(obj.stampSlice.loadCompVabc) = JVabc;

            % The frequency column comes from BOTH kernels. The linear one
            % owns every entry that is algebraic in the harmonic
            % coefficients; the nonlinear one owns the load torque's, which
            % can only be formed from the TIME-domain speed - writing it on
            % the coefficient instead is wrong by a power of N. Both halves
            % are state-dependent, so both reach the right-hand side; dW0b,
            % the bare additive w0 of the slip equation, still does not.
            dW0 = dW0 + reshape(dW0nlin,eqL,K1,nL);
            obj.stampVal(obj.stampSlice.loadCompW0) = dW0(:) + dW0b(:);

            % ---- Right-hand side ------------------------------------
            % Same accounting as STAMPIBRGFM; see the note there for why
            % dW0 reaches the right-hand side and dW0b does not.
            Jx = obj.blockJx(obj.map.loadCompNlin,r,JNlin);
            obj = obj.writeBlockRhs(r,eqL*K1, ...
                Jx - Fnlin + reshape(Flin,eqL*K1,nL) ...
                   + reshape(dW0,eqL*K1,nL)*w0);

            Vabc = obj.harmonicStackBus(obj.loadComp_bus,nL);
            JVb  = reshape(JVabc,4*K1,6*K1,nL);
            obj = obj.addBlockRhs(r, ...
                obj.harmonicLocalRows(idxHbLoadComp_Vs,eqL), ...
                reshape(pagemtimes(JVb,Vabc),4*K1,nL));
        end

        %% ------------- Harmonic gather/scatter helpers ----------------
        function x = devicePagesTd(obj,funGetIndices,nVar,nDev)
            % Inverse-DFT one DEVICE-LOCAL variable set for every device,
            % returned as [nVar x N x nDev] pages.
            %
            % The transpose before flattening matters: hbIdft splits a row
            % set into real and imaginary parts with (1:2:end)/(2:2:end),
            % so the rows it is handed must be INTERLEAVED (re,im) pairs.
            % base + locIdx comes out [nDev x nVarReIm], which flattens
            % device-fastest; transposing makes it variable-fastest, which
            % is the interleaving hbIdft needs and also puts the devices in
            % the order the reshape below assumes.
            v = obj.getTimeDomainSignal_iDFT(obj.tSamp, ...
                @(k)reshape(funGetIndices(k).',[],1));
            x = permute(reshape(v,nVar,nDev,obj.N),[1 3 2]);
        end

        function x = busPagesTd(obj,busList,nDev)
            % Inverse-DFT the three-phase voltage of each device's terminal
            % bus, as [3 x N x nDev] pages.
            %
            % getRowIdxAbc_fromBusNums is PHASE-major with the buses inner,
            % so hbIdft returns rows ordered (phase, bus) - bus fastest -
            % which is why the reshape is [nDev x 3] and not [3 x nDev].
            v = obj.getTimeDomainSignal_iDFT(obj.tSamp, ...
                @(k)obj.getRowIdxAbc_fromBusNums(busList,k));
            x = permute(reshape(v,nDev,3,obj.N),[2 3 1]);
        end

        function loc = harmonicLocalRows(obj,locRows,nEqns)
            % Local row indices of one variable across every harmonic,
            % harmonic-outer / variable-inner.
            loc = reshape((0:obj.K)*nEqns + locRows(:),[],1);
        end

        function S = harmonicStack(obj,rowBase,locRows,nEqns,nDev)
            % One variable's harmonic coefficients per device, as
            % [numel(locRows)*(K+1) x 1 x nDev] pages.
            loc = obj.harmonicLocalRows(locRows,nEqns);
            S = reshape(obj.x_v(rowBase(:).' + loc),numel(loc),1,nDev);
        end

        function S = harmonicStackBus(obj,busList,nDev)
            % The three-phase bus voltage across every harmonic, ordered to
            % match the kernel's dense block: harmonic-outer, node-inner.
            % Looped per device because getRowIdxAbc_fromBusNums is
            % phase-major across the WHOLE bus list, so one device's six
            % entries are not contiguous in it.
            S = zeros(6*(obj.K+1),1,nDev);
            for j = 1:nDev
                cols = reshape(obj.getRowIdxAbc_fromBusNum( ...
                    busList(j),0:obj.K),[],1);
                S(:,1,j) = obj.x_v(cols);
            end
        end

        function obj = writeHarmonicRows(obj,rowBase,locRows,nEqns,vals)
            % Overwrite one variable's rows across every harmonic. Used for
            % the rows a controller owns inside the machine's block.
            loc = obj.harmonicLocalRows(locRows,nEqns);
            rows = reshape(rowBase(:).' + loc,[],1);
            obj.I(rows) = vals(:);
        end
    end
end
