classdef EMT < handle
    properties (Constant)
        % System values
        f0 = 60
        w0 = 2*pi*EMT.f0 % Assume 60 Hz
        genrouP = 2 % Assume simple 2 pole generator

        % Equation setup
        numGenrouEqns = 17
        numExcDc1Eqns = 6 % 5 states + 1 vt calculation
        numTgov1dEqns = 3 % 2 states + 1 pmech calculation
        numIbrGFMEqns = 30  % 5 Park transf + 6 passive +
                           % 4 current control + 4 voltage controls +
                           % 5 droop + 6 virtual resistor
        numLoadCompEqns = 18 % IM:
                             % 4 stator (volt, curr) + 2 rotor (curr) +
                             % 2 speeds (wr,wslip) + 1 system theta +
                             % 3 iabc
                             % ZIP:
                             % 2 PQ + 2 Idq + 2 Vm (un-filt, filt)
                        % 2 Park V + 3 Park I

        absTolCheck = 1e-9

    end
    properties
        %% General properties

        % System
        baseMVA
        bus
        numBus = 0
        numGenrou = 0
        numExcDc1 = 0
        numTgov1d = 0
        numIbrGFM = 0
        numLoadComp = 0
        bus_BaseKv
        bus_IBase
        % Simulation (G*x_v = I, so I is negative)
        G
        I
        x_v
        x_v_prevTime
        tVec
        deltaT
        % Initial values
        x_v0
        x_v_prevTime0
        tx_I_prevTime0
        xfmr_I_prevTime0
        load_I_prevTime0
        % Indices for stamping
        % Triplet stamping (same contract as SALT): stampRow/stampCol are
        % built once by BUILDSTAMPMAPS, stampVal is refilled every Newton
        % iteration, and SALTASSEMBLEY turns the three into sparse G.
        stampRow
        stampCol
        stampVal
        stampSlice
        static_stampIdx
        static_stampJdx
        static_stampVals
        % Backward-Euler linear network values.  Shares static_stampIdx /
        % static_stampJdx with the trapezoidal set: stampLinear_G assembles
        % G_static for both methods and uses the STRUCTURAL UNION of the two
        % as the single sparsity pattern, so a coefficient that happens to
        % vanish under one method (e.g. the RL history term b when Rp=Inf)
        % cannot silently shift the pattern.
        %
        %
        static_numStampVals = 0
        genrou_numStampVals = 0
        excDc1_numStampVals = 0
        tgov1d_numStampVals = 0
        ibrGFM_numStampVals = 0
        loadComp_numStampVals = 0
        % Indices for constraining range
        idxThetaAll % Constrain to [0,2*pi)


        %% Static device properties

        % Txline - PI LINE (Br=branch,Sh=shunt)
        tx_busFrom
        tx_busTo
        tx_ReqBr
        tx_aHistBr
        tx_bHistBr
        tx_ReqSh
        tx_aHistSh
        tx_bHistSh
        % Backward-Euler companion coefficients.  Both sets are built at
        % setup so obj.intType can be flipped mid-simulation; the stamping
        % code picks the active set at call time (see getLinearCompanion).
        tx_I_prevTime
        tx_I_hist
        tx_R
        tx_L
        tx_C
        % Txmr - R-L
        xfmr_busFrom
        xfmr_busTo
        xfmr_Req
        xfmr_aHist
        xfmr_bHist
        xfmr_I_prevTime
        xfmr_I_hist
        xfmr_R
        xfmr_L
        % Const load - R-L
        load_bus
        load_Req
        load_aHist
        load_bHist
        load_I_prevTime
        load_I_hist
        load_R
        load_X
        load_Z

        %% Dynamic device properties
        % --------------------------- MACHINES ----------------------------

        % Generator (GENROU model)
        genrou_bus
        genrou_EsBase
        genrou_IsBase
        genrou_IsBaseRMS
        genrou_BaseMVA
        genrou_MW
        genrou_Mvar
        genrou_RMatrix
        genrou_LMatrix
        genrou_LspMatrix
        genrou_H
        genrou_KD
        genrou_TorqueM %**
        genrou_efd %**
        genrou_TorqueM0
        genrou_efd0
        %
        genrou_Lad
        genrou_Rfd
        genrou_Xq
        genrou_Xd
        genrou_Ra
        %

        % -------------------------- CONTROLLERS --------------------------
        % Controller parameters are per unit on the MACHINE base, not the
        % system base, and are used without conversion.

        % Exciter (EXDC1)
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
        excDc1_vref %**
        excDc1_vref0
        %
        excDc1_genToExc

        % Governor (TGOV1D)
        tgov1d_bus
        tgov1d_T1
        tgov1d_T2
        tgov1d_T3
        tgov1d_R
        tgov1d_Dt
        tgov1d_Pref %**
        tgov1d_Pref0
        %

        % --------------------------- IBRS --------------------------------
        % Droop-controlled inverter behind an LCL filter, with an inner
        % current loop, an outer voltage loop and a virtual resistor.
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
        ibrGFM_Mp
        ibrGFM_Mq
        ibrGFM_w_VR
        ibrGFM_Rv
        ibrGFM_pset
        ibrGFM_qset
        ibrGFM_vo_q_set

        % --------------------------- LOAD --------------------------------
        loadComp_bus
        loadComp_VBase
        loadComp_IBase
        loadComp_IsBaseRMS
        loadComp_MVA_base
        loadComp_MW
        loadComp_MVar
        loadCompIM_powerPercent
        loadCompIM_R
        loadCompIM_L
        loadCompIM_LTe
        loadCompIM_Tm0
        loadCompIM_m
        loadCompIM_H
        loadCompZip_PQabc
        loadCompZip_PQz
        loadCompZip_PQi
        loadCompZip_PQp
        loadCompZip_Tau

    end

    methods
        %% Class
        function obj = EMT(deltaT,t0,tEnd, baseMVA, bus_num)
            % Trapezoidal integration throughout - the machine, exciter and
            % governor companion models are only derived for Trap.

            % Time-domain model-index globals, used throughout the TD
            % stamping code. They belong to this class: a SALT solve that
            % never builds an EMT twin has no use for them, which is why
            % SALT's constructor populates only the steady-state set.
            emtInitIndices;

            obj.deltaT  = deltaT;
            obj.tVec    = t0:deltaT:tEnd;
            obj.baseMVA = baseMVA;
            obj.bus     = bus_num;
            obj.numBus  = numel(bus_num);
        end

        %% Simulation initializations

        function obj = systemValInitializations(obj, bus_Vm,bus_Va)
            % Sets x_v0 = x_v_prev = x_v = (vals from power flow)
            % Sets I_prev0 = I_prev = (calc from x_v_prev)

            % ******************** STATE VARIABLE *********************
            obj.x_v = zeros(obj.getSysSize,1);

            % Voltages
            % NOTE: Bus voltage of x_v0 is line-to-line (LL), RMS
            [bus_VComplex_3ph,bus_VComplex] = busValInitializations( ...
                bus_Vm,bus_Va);
            rowIdxBusAbc = obj.getRowIdxAbc_fromBusNums(obj.bus);
            assert(isequal(rowIdxBusAbc,(1:(3*obj.numBus)).'))
            obj.x_v(rowIdxBusAbc) = real(bus_VComplex_3ph);

            % Generators
            if obj.numGenrou > 0
                [obj,genrou_VAbsPhA,~,~,efd,Tm] = ...
                    genrouValInitializations(obj,bus_VComplex,bus_Va);
                obj.genrou_efd = efd;      obj.genrou_efd0 = efd;
                obj.genrou_TorqueM = Tm;   obj.genrou_TorqueM0 = Tm;
            end
            % Exciters
            if obj.numExcDc1 > 0
                [~,genrouIdx_excDc1] = ismembertol(obj.excDc1_bus,obj.genrou_bus);
                et = genrou_VAbsPhA(genrouIdx_excDc1); % Terminal voltage
                [obj,vref] = excDc1ValInitializations(obj,et);
                obj.excDc1_vref = vref;    obj.excDc1_vref0 = vref;
            end
            % Governors
            if obj.numTgov1d > 0
                [obj,pref] = tgov1dValInitializations(obj);
                obj.tgov1d_Pref = pref;    obj.tgov1d_Pref0 = pref;
            end
            % Grid forming IBRs
            if obj.numIbrGFM > 0
                [obj,pset,qset,vo_q_set] = ...
                    ibrGFMValInitializations(obj,bus_VComplex,bus_Va);
                obj.ibrGFM_pset = pset;
                obj.ibrGFM_qset = qset;
                obj.ibrGFM_vo_q_set = vo_q_set;
            end
            % Composite loads
            if obj.numLoadComp > 0
                [obj,Tm0,PQz,PQi,PQp] = loadCompValInitializations(obj, ...
                    bus_VComplex);
                obj.loadCompIM_Tm0 = Tm0;
                obj.loadCompZip_PQz = PQz;
                obj.loadCompZip_PQi = PQi;
                obj.loadCompZip_PQp = PQp;
            end

            % Set system frequency
            obj.x_v(obj.getRowIdx_wSys) = obj.w0;

            % Save initialization
            obj.x_v_prevTime = obj.x_v;
            obj.x_v0 = obj.x_v;
            obj.x_v_prevTime0 = obj.x_v_prevTime;

            % ******************** I_PREV VAR ********************
            % NOTE: This is NOT the same as "setIPrev_fromVPrev" because
            % initialization has to convert from phasor to real voltage and
            % current (setIPrev_fromVPrev uses history terms, which do not
            % exist yet at initialization, so the reactances are used).
            [tx_I_prevTime_complex,xfmr_I_prevTime_complex, ...
                load_I_prevTime_complex] = ...
                obj.historyCurrentInitializations(bus_VComplex_3ph,obj.w0);
            obj.tx_I_prevTime   = real(tx_I_prevTime_complex);
            obj.xfmr_I_prevTime = real(xfmr_I_prevTime_complex);
            obj.load_I_prevTime = real(load_I_prevTime_complex);

            obj.tx_I_prevTime0   = obj.tx_I_prevTime;
            obj.xfmr_I_prevTime0 = obj.xfmr_I_prevTime;
            obj.load_I_prevTime0 = obj.load_I_prevTime;
        end

        function [tx_I_prev_complex,xfmr_I_prev_complex,load_I_prev_complex] = ...
                historyCurrentInitializations(obj,bus_VComplex_3ph,w)

            % Set I_prev for Txline
            tx_VfComplex = bus_VComplex_3ph( ...
                obj.getRowIdxAbc_fromBusNums(obj.tx_busFrom));
            tx_VtComplex = bus_VComplex_3ph( ...
                obj.getRowIdxAbc_fromBusNums(obj.tx_busTo));
            Ibranch_complex = (tx_VfComplex - tx_VtComplex)./ ...
                repmat(obj.tx_R + 1j*w*obj.tx_L, 3,1);
            IShFrom_complex = tx_VfComplex.*repmat(1j*w*obj.tx_C, 3,1);
            IShTo_complex   = tx_VtComplex.*repmat(1j*w*obj.tx_C, 3,1);
            tx_I_prev_complex = [Ibranch_complex,IShFrom_complex,IShTo_complex];
            % Verify shape: 2 shunt currents + 1 branch current
            if ~isempty(obj.tx_busFrom)
                assert(isequal(size(tx_I_prev_complex),[3*numel(obj.tx_busFrom),3]));
            end

            % Set I_prev for transformer
            xfmr_VfComplex = bus_VComplex_3ph( ...
                obj.getRowIdxAbc_fromBusNums(obj.xfmr_busFrom));
            xfmr_VtComplex = bus_VComplex_3ph( ...
                obj.getRowIdxAbc_fromBusNums(obj.xfmr_busTo));
            xfmr_I_prev_complex = (xfmr_VfComplex - xfmr_VtComplex)./ ...
                repmat(obj.xfmr_R + 1j*w*obj.xfmr_L, 3,1);
            if ~isempty(obj.xfmr_busFrom)
                assert(length(xfmr_I_prev_complex)==3*numel(obj.xfmr_busFrom))
            end

            % Set I_prev for load
            if isempty(obj.load_bus)
                load_I_prev_complex = [];
            elseif w==obj.w0
                load_I_prev_complex = bus_VComplex_3ph( ...
                    obj.getRowIdxAbc_fromBusNums(obj.load_bus) ) ...
                    ./ repmat(obj.load_Z,3,1);
            else
                load_I_prev_complex = zeros(3*numel(obj.load_bus),1);
                % Model loads as RL?
                idxRLSeries = obj.load_X > 0;
                if any(idxRLSeries,'all')
                    V = bus_VComplex_3ph( obj.getRowIdxAbc_fromBusNums( ...
                        obj.load_bus(idxRLSeries)) );
                    R = obj.load_R(idxRLSeries);
                    L = obj.load_X(idxRLSeries) / obj.w0;
                    load_I_prev_complex(repmat(idxRLSeries,3,1)) = ...
                        V .* repmat(1./(R + 1j*w*L), 3,1);
                end
                % Model as R//C?
                idxRCPar = obj.load_X < 0;
                if any(idxRCPar)
                    V = bus_VComplex_3ph( obj.getRowIdxAbc_fromBusNums( ...
                        obj.load_bus(idxRCPar)) );
                    R = 1 ./ real(1./obj.load_Z(idxRCPar));
                    C = imag(1./obj.load_Z(idxRCPar)) / obj.w0;
                    load_I_prev_complex(repmat(idxRCPar,3,1)) = ...
                        V .* repmat(1./R + 1j*w*C, 3,1);
                end
                % Model as R?
                idxR = ~(idxRLSeries+idxRCPar);
                if any(idxR)
                    V = bus_VComplex_3ph( obj.getRowIdxAbc_fromBusNums( ...
                        obj.load_bus(idxR)) );
                    load_I_prev_complex(repmat(idxR,3,1)) = ...
                        V .* repmat(1./obj.load_R(idxR), 3,1);
                end
            end
        end

        function boundThetaRadiansInitializations(obj)
            % Collect every rotating-frame angle so it can be wrapped into
            % [0,2*pi) after each time step, which keeps the Park
            % transforms well conditioned over long runs.
            idxGenrouTheta = []; idxIbrGFMTheta = []; idxLoadCompTheta = [];
            if obj.numGenrou>0
                [~,~,~,~,genrouIdxTheta] = getGenrouIndicesTD;
                idxGenrouTheta = obj.getRowIdx_fromGenrouIndices(1:obj.numGenrou)-1 ...
                    + genrouIdxTheta;
            end
            if obj.numIbrGFM>0
                [~,~,~,~,~,~,~,~,~,~,~,~,ibrGFMIdxTheta] = getIbrGFMIndicesTD;
                idxIbrGFMTheta = obj.getRowIdx_fromIbrGFMIndices(1:obj.numIbrGFM) ...
                    -1 + ibrGFMIdxTheta;
            end
            if obj.numLoadComp>0
                [~,~,~,idxThetaSys] = getLoadCompIndicesTD;
                idxLoadCompTheta = ...
                    obj.getRowIdx_fromLoadCompIndices(1:obj.numLoadComp)-1 ...
                    + idxThetaSys;
            end
            obj.idxThetaAll = [idxGenrouTheta,idxIbrGFMTheta,idxLoadCompTheta];
        end

        %% Simulation stamping and upating

        function obj = buildStampMaps(obj)
            %BUILDSTAMPMAPS  Precompute the (row,col) of every triplet.
            %
            % Same contract as SALT.buildStampMaps: obj.stampRow and
            % obj.stampCol are fixed for the whole run and only
            % obj.stampVal is refilled each Newton iteration. The linear
            % network contributes constant companion conductances, so its
            % values are computed once here and cached in
            % obj.static_stampVals.

            % ---------------- Linear network -----------------------------
            Rs = {}; Cs = {}; Vs = {};
            function addBlock(rowIdx,colIdx,block)
                [rr,cc] = saltStampTriplets(rowIdx,colIdx);
                Rs{end+1} = rr; Cs{end+1} = cc; Vs{end+1} = block(:); %#ok<AGROW>
            end

            % Txline: RL series branch plus a capacitive shunt at each end
            for idx = 1:numel(obj.tx_busFrom)
                idxF = obj.getRowIdxAbc_fromBusNum(obj.tx_busFrom(idx));
                idxT = obj.getRowIdxAbc_fromBusNum(obj.tx_busTo(idx));
                idxFT = reshape([idxF,idxT].',[],1);
                addBlock(idxFT,idxFT, ...
                    obj.stampOfTxGMat(obj.tx_ReqBr(idx),obj.tx_ReqSh(idx)));
            end
            % XFMR: RL series branch
            for idx = 1:numel(obj.xfmr_busFrom)
                idxF = obj.getRowIdxAbc_fromBusNum(obj.xfmr_busFrom(idx));
                idxT = obj.getRowIdxAbc_fromBusNum(obj.xfmr_busTo(idx));
                idxFT = reshape([idxF,idxT].',[],1);
                addBlock(idxFT,idxFT,obj.stampOfXfmrGMat(obj.xfmr_Req(idx)));
            end
            % Constant-impedance load
            for idx = 1:numel(obj.load_bus)
                idxF = obj.getRowIdxAbc_fromBusNum(obj.load_bus(idx));
                addBlock(idxF,idxF,obj.stampOfLoadGMat(obj.load_Req(idx)));
            end

            % Current injected by the nonlinear devices.
            % Current direction: machines and inverters produce power, so
            % their current flows out of the device into the circuit and is
            % stamped negative. Loads consume power, so their current flows
            % into the device and is stamped positive.
            [~,~,~,idxGenrouWr,~,idxIabc] = getGenrouIndicesTD;
            for idx = 1:obj.numGenrou
                addBlock(obj.getRowIdxAbc_fromBusNum(obj.genrou_bus(idx)), ...
                    obj.getRowIdx_fromGenrouIndices(idx)-1+idxIabc,-eye(3));
            end
            [idxIabcIbrGFM,~,~,~,~,~,~,~,~,idxIbrGFMWr] = getIbrGFMIndicesTD;
            for idx = 1:obj.numIbrGFM
                addBlock(obj.getRowIdxAbc_fromBusNum(obj.ibrGFM_bus(idx)), ...
                    obj.getRowIdx_fromIbrGFMIndices(idx)-1+idxIabcIbrGFM,-eye(3));
            end
            if obj.numLoadComp>0
                [~,idxIabcLd] = getLoadCompIndicesTD;
                for idx = 1:obj.numLoadComp
                    addBlock(obj.getRowIdxAbc_fromBusNum( ...
                        obj.loadComp_bus(idx)), ...
                        obj.getRowIdx_fromLoadCompIndices(idx)-1+idxIabcLd, ...
                        eye(3));
                end
            end

            % System frequency: an MVA-weighted average of every rotating
            % machine and inverter speed.
            idxWSys = obj.getRowIdx_wSys;
            idxGenrou_wr = idxGenrouWr + ...
                obj.getRowIdx_fromGenrouIndices(1:obj.numGenrou)-1;
            idxIbrGFM_ws = idxIbrGFMWr + ...
                obj.getRowIdx_fromIbrGFMIndices(1:obj.numIbrGFM)-1;
            addBlock(idxWSys,[idxGenrou_wr,idxIbrGFM_ws,idxWSys], ...
                obj.stampOfWsysGMat);

            obj.static_stampIdx  = vertcat(Rs{:});
            obj.static_stampJdx  = vertcat(Cs{:});
            obj.static_stampVals = vertcat(Vs{:});
            obj.static_numStampVals = numel(obj.static_stampIdx);

            % ---------------- Nonlinear device blocks --------------------
            % Each model contributes a fixed per-device pattern; the flat
            % stamping routines fill the matching value slots.
            if obj.numGenrou > 0
                [idx,jdx,idx_dFParkVdq_dvabc,jdx_dFParkVdq_dvabc] = ...
                    calcGenrouIndices;
                idxGenrouRow = obj.getRowIdx_fromGenrouIndices(1:obj.numGenrou)-1;
                idxVabcRow_hasGenrou = obj.getRowIdxAbc_fromBusNums(obj.genrou_bus.');
                genrou_stampIdx = [idxGenrouRow+idx; idxGenrouRow+idx_dFParkVdq_dvabc];
                genrou_stampJdx = [idxGenrouRow+jdx; ...
                    idxVabcRow_hasGenrou(1,:)+jdx_dFParkVdq_dvabc(1:3); ...
                    idxVabcRow_hasGenrou(2,:)+jdx_dFParkVdq_dvabc(4:6); ...
                    idxVabcRow_hasGenrou(3,:)+jdx_dFParkVdq_dvabc(7:9)];
                assert(numel(genrou_stampIdx)==numel(genrou_stampJdx))
                obj.genrou_numStampVals = numel(genrou_stampIdx)/obj.numGenrou;
                % Every machine must have both an exciter and a governor
                assert(isequal(sort(obj.genrou_bus),sort(obj.excDc1_bus)))
                assert(isequal(sort(obj.genrou_bus),sort(obj.tgov1d_bus)))
            else
                genrou_stampIdx = []; genrou_stampJdx = [];
            end
            if obj.numExcDc1 > 0
                [idx,jdx, ...
                    idx_dFefd_dx1,jdx_dFefd_dx1, ...
                    idx_dFefd_defd,jdx_dFefd_defd, ...
                    idx_dFefd_dwr,jdx_dFefd_dwr, ...
                    idx_dFexcVt_dvt,jdx_dFexcVt_dvt, ...
                    idx_dFexcVt_dvabc,jdx_dFexcVt_dvabc] = calcExcDc1Indices;
                idxExcDc1Row = obj.getRowIdx_fromExcDc1Indices(1:obj.numExcDc1)-1;
                [~,genrouIdx_hasExcDc1] = ismembertol(obj.genrou_bus,obj.excDc1_bus);
                idxGenrouRow_hasExcDc1 = idxGenrouRow(genrouIdx_hasExcDc1);
                idxVabcRow_hasExcDc1 = obj.getRowIdxAbc_fromBusNums(obj.excDc1_bus.');
                excDc1_stampIdx = [idxExcDc1Row+idx; ...
                    idxGenrouRow_hasExcDc1+idx_dFefd_dx1; ...
                    idxGenrouRow_hasExcDc1+idx_dFefd_defd; ...
                    idxGenrouRow_hasExcDc1+idx_dFefd_dwr; ...
                    idxExcDc1Row+idx_dFexcVt_dvt; ...
                    idxExcDc1Row+idx_dFexcVt_dvabc];
                excDc1_stampJdx = [idxExcDc1Row+jdx; ...
                    idxExcDc1Row+jdx_dFefd_dx1; ...
                    idxGenrouRow_hasExcDc1+jdx_dFefd_defd; ...
                    idxGenrouRow_hasExcDc1+jdx_dFefd_dwr; ...
                    idxExcDc1Row+jdx_dFexcVt_dvt; ...
                    idxVabcRow_hasExcDc1(1,:)+jdx_dFexcVt_dvabc(1); ...
                    idxVabcRow_hasExcDc1(2,:)+jdx_dFexcVt_dvabc(2); ...
                    idxVabcRow_hasExcDc1(3,:)+jdx_dFexcVt_dvabc(3)];
                assert(numel(excDc1_stampIdx)==numel(excDc1_stampJdx))
                obj.excDc1_numStampVals = numel(excDc1_stampIdx)/obj.numExcDc1;
            else
                excDc1_stampIdx = []; excDc1_stampJdx = [];
            end
            if obj.numTgov1d > 0
                [idx,jdx, ...
                    idx_dFx2_dwr,jdx_dFx2_dwr, ...
                    idx_dFx3_dwr,jdx_dFx3_dwr, ...
                    idx_dFTm_dPmech,jdx_dFTm_dPmech, ...
                    idx_dFTm_dwr,jdx_dFTm_dwr, ...
                    idx_dFTm_dTm,jdx_dFTm_dTm] = calcTgov1dIndices;
                idxTgov1dRow = obj.getRowIdx_fromTgov1dIndices(1:obj.numTgov1d)-1;
                [~,genrouIdx_hasTgov1d] = ismembertol(obj.genrou_bus,obj.tgov1d_bus);
                idxGenrouRow_hasTgov1d = idxGenrouRow(genrouIdx_hasTgov1d);
                tgov1d_stampIdx = [idxTgov1dRow+idx; ...
                    idxTgov1dRow+idx_dFx2_dwr; idxTgov1dRow+idx_dFx3_dwr; ...
                    idxGenrouRow_hasTgov1d+idx_dFTm_dPmech; ...
                    idxGenrouRow_hasTgov1d+idx_dFTm_dwr; ...
                    idxGenrouRow_hasTgov1d+idx_dFTm_dTm];
                tgov1d_stampJdx = [idxTgov1dRow+jdx; ...
                    idxGenrouRow+jdx_dFx2_dwr; idxGenrouRow+jdx_dFx3_dwr; ...
                    idxTgov1dRow+jdx_dFTm_dPmech; ...
                    idxGenrouRow_hasTgov1d+jdx_dFTm_dwr; ...
                    idxGenrouRow_hasTgov1d+jdx_dFTm_dTm];
                assert(numel(tgov1d_stampIdx)==numel(tgov1d_stampJdx))
                obj.tgov1d_numStampVals = numel(tgov1d_stampIdx)/obj.numTgov1d;
            else
                tgov1d_stampIdx = []; tgov1d_stampJdx = [];
            end
            if obj.numIbrGFM > 0
                [idx,jdx,idx_dFParkVdq_dvabc,jdx_dFParkVdq_dvabc] = ...
                    calcIbrGFMIndices;
                idxIbrGFMRow = obj.getRowIdx_fromIbrGFMIndices(1:obj.numIbrGFM)-1;
                idxVabcRow_hasIbrGFM = obj.getRowIdxAbc_fromBusNums( ...
                    obj.ibrGFM_bus(:).');
                ibrGFM_stampIdx = [idxIbrGFMRow+idx; ...
                    idxIbrGFMRow+idx_dFParkVdq_dvabc];
                ibrGFM_stampJdx = [idxIbrGFMRow+jdx; ...
                    idxVabcRow_hasIbrGFM(1,:)+jdx_dFParkVdq_dvabc(1:2); ...
                    idxVabcRow_hasIbrGFM(2,:)+jdx_dFParkVdq_dvabc(3:4); ...
                    idxVabcRow_hasIbrGFM(3,:)+jdx_dFParkVdq_dvabc(5:6)];
                assert(numel(ibrGFM_stampIdx)==numel(ibrGFM_stampJdx))
                obj.ibrGFM_numStampVals = numel(ibrGFM_stampIdx)/obj.numIbrGFM;
            else
                ibrGFM_stampIdx = []; ibrGFM_stampJdx = [];
            end
            if obj.numLoadComp > 0
                [idx,jdx, ...
                    idx_dFParkVdq_dvabc,jdx_dFParkVdq_dvabc, ...
                    idx_dFslip_dwSys,jdx_dFslip_dwSys, ...
                    idx_dFsys_dwSys,jdx_dFsys_dwSys, ...
                    idx_dFs_dwSys,jdx_dFs_dwSys, ...
                    idx_dFswing_dwSys,jdx_dFswing_dwSys] = ...
                    calcLoadCompIndices;
                idxLoadCompRow = obj.getRowIdx_fromLoadCompIndices( ...
                    1:obj.numLoadComp)-1;
                idxVabcRow_hasCompRow = obj.getRowIdxAbc_fromBusNums( ...
                    obj.loadComp_bus.');
                idxWSysRow = obj.getRowIdx_wSys*ones(1,obj.numLoadComp);
                loadComp_stampIdx = [idxLoadCompRow+idx; ...
                    idxLoadCompRow+idx_dFParkVdq_dvabc; ...
                    idxLoadCompRow+idx_dFslip_dwSys; ...
                    idxLoadCompRow+idx_dFsys_dwSys; ...
                    idxLoadCompRow+idx_dFs_dwSys; ...
                    idxLoadCompRow+idx_dFswing_dwSys];
                loadComp_stampJdx = [idxLoadCompRow+jdx; ...
                    idxVabcRow_hasCompRow(1,:)+jdx_dFParkVdq_dvabc(1:2); ...
                    idxVabcRow_hasCompRow(2,:)+jdx_dFParkVdq_dvabc(3:4); ...
                    idxVabcRow_hasCompRow(3,:)+jdx_dFParkVdq_dvabc(5:6); ...
                    idxWSysRow+jdx_dFslip_dwSys; ...
                    idxWSysRow+jdx_dFsys_dwSys; ...
                    idxWSysRow+jdx_dFs_dwSys; ...
                    idxWSysRow+jdx_dFswing_dwSys];
                assert(numel(loadComp_stampIdx)==numel(loadComp_stampJdx))
                obj.loadComp_numStampVals = numel(loadComp_stampIdx)/ ...
                    obj.numLoadComp;
            else
                loadComp_stampIdx = []; loadComp_stampJdx = [];
            end

            % ---------------- Scatter into the triplet arrays ------------
            numStampVals = obj.getStampSize;
            obj.stampRow = zeros(numStampVals,1);
            obj.stampCol = zeros(numStampVals,1);
            obj.stampVal = zeros(numStampVals,1);

            obj.stampRow(obj.getStampIdx_static) = obj.static_stampIdx;
            obj.stampCol(obj.getStampIdx_static) = obj.static_stampJdx;
            obj.stampVal(obj.getStampIdx_static) = obj.static_stampVals;

            slot = obj.getStampIdx_fromGenrouIndices(1:obj.numGenrou)-1 ...
                + (1:obj.genrou_numStampVals).';
            obj.stampRow(slot) = genrou_stampIdx;
            obj.stampCol(slot) = genrou_stampJdx;
            obj.stampSlice.genrou = slot(:);

            slot = obj.getStampIdx_fromExcDc1Indices(1:obj.numExcDc1)-1 ...
                + (1:obj.excDc1_numStampVals).';
            obj.stampRow(slot) = excDc1_stampIdx;
            obj.stampCol(slot) = excDc1_stampJdx;
            obj.stampSlice.excDc1 = slot(:);

            slot = obj.getStampIdx_fromTgov1dIndices(1:obj.numTgov1d)-1 ...
                + (1:obj.tgov1d_numStampVals).';
            obj.stampRow(slot) = tgov1d_stampIdx;
            obj.stampCol(slot) = tgov1d_stampJdx;
            obj.stampSlice.tgov1d = slot(:);

            slot = obj.getStampIdx_fromIbrGFMIndices(1:obj.numIbrGFM)-1 ...
                + (1:obj.ibrGFM_numStampVals).';
            obj.stampRow(slot) = ibrGFM_stampIdx;
            obj.stampCol(slot) = ibrGFM_stampJdx;
            obj.stampSlice.ibrGFM = slot(:);

            slot = obj.getStampIdx_fromLoadCompIndices(1:obj.numLoadComp)-1 ...
                + (1:obj.loadComp_numStampVals).';
            obj.stampRow(slot) = loadComp_stampIdx;
            obj.stampCol(slot) = loadComp_stampJdx;
            obj.stampSlice.loadComp = slot(:);

            idxStatic = obj.getStampIdx_static;
            obj.stampSlice.static = idxStatic(:);
        end


        function GSysFreq = stampOfWsysGMat(obj)
            % F = SUM{ baseMVA[i]*f[i] } / SUM{ baseMVA[i] }

            % Indices in order idxGenrou_wr, idxIbrGFM_ws, idxWSys
            baseVals = [obj.genrou_BaseMVA;obj.ibrGFM_MVA_base].';
            sumBase = sum(baseVals);
            GSysFreq = [baseVals/sumBase,-1];


        end


        function obj = stampNonlinear_G_I(obj)
            % Refill the triplet VALUES for the current Newton guess and
            % rebuild G. obj.stampRow / obj.stampCol were fixed once by
            % BUILDSTAMPMAPS, exactly as in SALT.
            all_stampVal = zeros(numel(obj.stampRow),1);
            all_stampVal(obj.getStampIdx_static) = obj.static_stampVals;

            % Local copies of the state and of the RHS. Every vectorized
            % block below reads x_v/x_p and writes I_v; obj.I is written
            % back once, at the end. This keeps the property access out of
            % the inner code.
            x_v = obj.x_v;
            x_p = obj.x_v_prevTime;
            I_v = obj.I;

            % The machine, exciter, governor, grid-forming IBR and
            % composite-load blocks are flat, fully vectorized routines:
            % they gather the whole device population from x_v at once and
            % write straight into all_stampVal / I_v - no per-device loop.
            % They are FUNCTIONS rather than scripts on purpose: a script
            % runs in a dynamic workspace that MATLAB cannot JIT-compile,
            % which cost ~4-6x here and swamped any saving from avoiding
            % the call.
            % ORDER MATTERS: the exciter overwrites the machine's Efd row
            % and the governor overwrites its Tm row, exactly as the
            % original per-device loop ordering did.
            if obj.numGenrou > 0
                [all_stampVal,I_v] = stampOfGenrouMatrices_TD_flat( ...
                    obj,all_stampVal,I_v,x_v,x_p);
            end
            if obj.numExcDc1 > 0
                [all_stampVal,I_v] = stampOfExcDc1Matrices_TD_flat( ...
                    obj,all_stampVal,I_v,x_v,x_p);
            end
            if obj.numTgov1d > 0
                [all_stampVal,I_v] = stampOfTgov1dMatrices_TD_flat( ...
                    obj,all_stampVal,I_v,x_v,x_p);
            end
            if obj.numIbrGFM > 0
                [all_stampVal,I_v] = stampOfIbrGFMMatrices_TD_flat( ...
                    obj,all_stampVal,I_v,x_v,x_p);
            end
            if obj.numLoadComp > 0
                [all_stampVal,I_v] = ...
                    stampOfCompositeLoadMatrices_TD_flat( ...
                    obj,all_stampVal,I_v,x_v,x_p);
            end

            % Write back and assemble from the triplets
            obj.I = I_v;
            obj.stampVal = all_stampVal;
            obj.G = saltAssembleY(obj.stampRow,obj.stampCol, ...
                obj.stampVal,obj.getSysSize);
        end

        function obj = stampLinear_I(obj,idxTvec) %#ok<INUSD>
            % Represents stamping history current -- is the same in each NR
            % iteration
            % NOTE: G*x_v = I, so I is negative (i.e. represents current
            %       going into the node)

            sysSize = obj.getSysSize;
            obj.I = zeros(sysSize,1);

            % Companion coefficients for the ACTIVE integration method
            comp = obj.getLinearCompanion;

            % Stamp Tx
            numTx = numel(obj.tx_busFrom);
            idxF = obj.getRowIdxAbc_fromBusNums(obj.tx_busFrom.');
            idxT = obj.getRowIdxAbc_fromBusNums(obj.tx_busTo.');
            for idx = 1:numel(obj.tx_busFrom)
                idxDevice = idx + [0;numTx;2*numTx];
                % Stamp
                [IHistBr,IHistShFrom,IHistShTo] = obj.stampOfTxIhist( ...
                    obj.tx_I_prevTime(idxDevice,1), ...
                    obj.tx_I_prevTime(idxDevice,2), ...
                    obj.tx_I_prevTime(idxDevice,3), ...
                    obj.x_v_prevTime(idxF(:,idx)), ...
                    obj.x_v_prevTime(idxT(:,idx)), ...
                    comp.txAHistBr(idx),comp.txBHistBr(idx), ...
                    comp.txAHistSh(idx),comp.txBHistSh(idx));
                obj.I(idxF(:,idx)) = obj.I(idxF(:,idx)) - (IHistBr + IHistShFrom);
                obj.I(idxT(:,idx)) = obj.I(idxT(:,idx)) - (-IHistBr + IHistShTo);
                % History
                obj.tx_I_hist(idxDevice,:) = [IHistBr,IHistShFrom,IHistShTo];
            end

            % Stamp XFMR
            numXfmr = numel(obj.xfmr_busFrom);
            idxF = obj.getRowIdxAbc_fromBusNums(obj.xfmr_busFrom.');
            idxT = obj.getRowIdxAbc_fromBusNums(obj.xfmr_busTo.');
            for idx = 1:numXfmr
                idxDevice = idx + [0;numXfmr;2*numXfmr];
                % Stamp
                IHist = obj.stampOfXfmrIhist(obj.xfmr_I_prevTime(idxDevice), ...
                    obj.x_v_prevTime(idxF(:,idx)),obj.x_v_prevTime(idxT(:,idx)), ...
                    comp.xfmrAHist(idx),comp.xfmrBHist(idx));
                obj.I(idxF(:,idx)) = obj.I(idxF(:,idx)) - (IHist);
                obj.I(idxT(:,idx)) = obj.I(idxT(:,idx)) - (-IHist);
                % History
                obj.xfmr_I_hist(idxDevice) = IHist;
            end

            % Stamp Load
            numLoad = numel(obj.load_bus);
            idxF = obj.getRowIdxAbc_fromBusNums(obj.load_bus.');
            for idx = 1:numLoad
                idxDevice = idx + [0;numLoad;2*numLoad];
                % Stamp
                IHist = obj.stampOfLoadIhist(obj.load_I_prevTime(idxDevice), ...
                    obj.x_v_prevTime(idxF(:,idx)), ...
                    comp.loadAHist(idx),comp.loadBHist(idx));
                obj.I(idxF(:,idx)) = obj.I(idxF(:,idx)) - (IHist);
                % History
                obj.load_I_hist(idxDevice) = IHist;
            end
            % Stamp C at each node

            % Stamp higher harmonics as current source
        end

        function obj = setIPrev_fromVPrev(obj)

            % Companion coefficients for the ACTIVE integration method
            comp = obj.getLinearCompanion;

            % Set I_prev for Txline
            tx_Vf = obj.x_v_prevTime(obj.getRowIdxAbc_fromBusNums(obj.tx_busFrom));
            tx_Vt = obj.x_v_prevTime(obj.getRowIdxAbc_fromBusNums(obj.tx_busTo));
            Ibranch = (tx_Vf - tx_Vt)./repmat(comp.txReqBr,3,1) ...
                + obj.tx_I_hist(:,1);
            IShFrom = tx_Vf./repmat(comp.txReqSh,3,1) ...
                + obj.tx_I_hist(:,2);
            IShTo = tx_Vt./repmat(comp.txReqSh,3,1) ...
                + obj.tx_I_hist(:,3);
            % Set
            obj.tx_I_prevTime = [Ibranch,IShFrom,IShTo];

            % Set I_prev for transformer
            xfmr_Vf = obj.x_v_prevTime(obj.getRowIdxAbc_fromBusNums(obj.xfmr_busFrom));
            xfmr_Vt = obj.x_v_prevTime(obj.getRowIdxAbc_fromBusNums(obj.xfmr_busTo));
            xfmrIbranch = (xfmr_Vf - xfmr_Vt)./repmat(comp.xfmrReq,3,1) ...
                + obj.xfmr_I_hist;
            % Set
            obj.xfmr_I_prevTime = xfmrIbranch;

            % Set I_prev for load
            obj.load_I_prevTime = ....
                obj.x_v_prevTime( obj.getRowIdxAbc_fromBusNums(obj.load_bus) ) ...
                ./ repmat(comp.loadReq,3,1) ...
                + obj.load_I_hist;

            % Set I_prev for C at every node
        end

        function c = getLinearCompanion(obj)
            % Trapezoidal companion coefficients for the linear network.
            c.txReqBr  = obj.tx_ReqBr;   c.txAHistBr = obj.tx_aHistBr;
            c.txBHistBr= obj.tx_bHistBr; c.txReqSh   = obj.tx_ReqSh;
            c.txAHistSh= obj.tx_aHistSh; c.txBHistSh = obj.tx_bHistSh;
            c.xfmrReq  = obj.xfmr_Req;   c.xfmrAHist = obj.xfmr_aHist;
            c.xfmrBHist= obj.xfmr_bHist;
            c.loadReq  = obj.load_Req;   c.loadAHist = obj.load_aHist;
            c.loadBHist= obj.load_bHist;
        end

        function boundThetaRadians(obj)
            obj.x_v(obj.idxThetaAll) = mod(obj.x_v(obj.idxThetaAll),2*pi);

        end


        %% Simulation helper functions
        %% Set device params

        % ***** TXLINE PARAMS *****
        function obj = setTxParams(obj,tx_busFrom,tx_busTo, tx_R,tx_L,tx_C)
            % 3ph PI line: an RL series branch with a capacitive shunt at
            % each end, each replaced by its trapezoidal companion model.
            obj.tx_busFrom = tx_busFrom;
            obj.tx_busTo = tx_busTo;
            obj.tx_I_hist = zeros(3*numel(obj.tx_busFrom),3);
            obj.tx_R = tx_R;
            obj.tx_L = tx_L;
            obj.tx_C = tx_C;
            % Damping resistors keep the companion models well conditioned
            tx_RParallel_L = obj.getDampingR_ParL(tx_L,obj.deltaT);
            tx_RSeries_C   = obj.getDampingR_SeriesC(tx_C,obj.deltaT);
            [obj.tx_ReqBr,obj.tx_aHistBr,obj.tx_bHistBr] = ...
                obj.calcReqHistParams_RLSeriesWithDampingRp(...
                obj.tx_R,obj.tx_L,tx_RParallel_L,obj.deltaT);
            [obj.tx_ReqSh,obj.tx_aHistSh,obj.tx_bHistSh] = ....
                obj.calcReqHistParams_CWithDampingRs(...
                obj.tx_C,tx_RSeries_C,obj.deltaT);
        end

         % ***** XFMR PARAMS *****
        function obj = setXfmrParams(obj,xfmr_busFrom,xfmr_busTo,xfmr_R,xfmr_L)
            % 3ph R-L series branch, trapezoidal companion model.
            obj.xfmr_busFrom = xfmr_busFrom;
            obj.xfmr_busTo = xfmr_busTo;
            obj.xfmr_I_hist = zeros(3*numel(obj.xfmr_busFrom),1);
            obj.xfmr_R = xfmr_R;
            obj.xfmr_L = xfmr_L;
            % Damping resistor
            xfmr_RParallel_L = obj.getDampingR_ParL(xfmr_L,obj.deltaT);
            [obj.xfmr_Req,obj.xfmr_aHist,obj.xfmr_bHist] = ...
                obj.calcReqHistParams_RLSeriesWithDampingRp(...
                obj.xfmr_R,obj.xfmr_L,xfmr_RParallel_L,obj.deltaT);
        end

         % ***** CONST LOAD PARAMS *****
        function obj = setLoadParams(obj,load_bus, load_Z)
            % Constant-impedance load, modelled as R-L, R//C or plain R
            % depending on the sign of the reactance.
            obj.load_bus = load_bus;
            numLoadBus = numel(load_bus);
            obj.load_I_hist = zeros(3*numLoadBus,1);
            obj.load_Z = load_Z;
            obj.load_X = imag(obj.load_Z);
            obj.load_R = real(obj.load_Z);
            obj.load_Req   = zeros(numLoadBus,1);
            obj.load_aHist = zeros(numLoadBus,1);
            obj.load_bHist = zeros(numLoadBus,1);
            % Model loads as RL?
            idxRLSeries = obj.load_X > 0;
            if any(idxRLSeries,'all')
                R = obj.load_R(idxRLSeries);
                L = obj.load_X(idxRLSeries) / obj.w0;
                RParallel_L = obj.getDampingR_ParL(L,obj.deltaT);
                [obj.load_Req(idxRLSeries),obj.load_aHist(idxRLSeries), ...
                    obj.load_bHist(idxRLSeries)] = ...
                    obj.calcReqHistParams_RLSeriesWithDampingRp(...
                    R,L,RParallel_L,obj.deltaT);
            end
            % Model as R//C?
            idxRCPar = obj.load_X < 0;
            if any(idxRCPar)
                R = 1 ./ real(1./load_Z(idxRCPar));
                C = imag(1./obj.load_Z(idxRCPar)) / obj.w0;
                RSeries_C = obj.getDampingR_SeriesC(C,obj.deltaT);
                [obj.load_Req(idxRCPar),obj.load_aHist(idxRCPar), ...
                    obj.load_bHist(idxRCPar)] = ...
                    EMT.calcReqHistParams_RParCWithDampingRs( ...
                    R,C,RSeries_C,obj.deltaT);
            end
            % Model as R?
            idxR = ~(idxRLSeries+idxRCPar);
            if any(idxR)
                obj.load_Req(idxR) = obj.load_R(idxR);
                obj.load_aHist(idxR) = 0;
                obj.load_bHist(idxR) = 0;
            end
        end

         % ***** COMPOSITE LOAD PARAMS *****
         function obj = setLoadCompParams(obj,loadComp_bus, ...
                 ... % Base
                 loadComp_VBase,loadComp_IBase,loadComp_MVA_base, ...
                 ... % PQ
                 loadComp_MW,loadComp_MVar, ...
                 ... % Parameters
                 loadCompIM_powerPercent, ...
                 loadCompIM_Rs,loadCompIM_Rr, ...
                 loadCompIM_Lss,loadCompIM_Lrr,loadCompIM_Lm, ...
                 loadCompIM_H, ...
                 loadCompIM_m, ...
                 loadCompZip_PQabc, ...
                 loadCompZip_Tau)
             % Set bus numbers
             obj.loadComp_bus = loadComp_bus;
             obj.numLoadComp = numel(loadComp_bus);
             % Set base
             obj.loadComp_VBase = loadComp_VBase;
             obj.loadComp_IBase = loadComp_IBase;
             obj.loadComp_IsBaseRMS = loadComp_IBase/sqrt(2);
             obj.loadComp_MVA_base = loadComp_MVA_base;
             % Set power
             obj.loadComp_MW = loadComp_MW;
             obj.loadComp_MVar = loadComp_MVar;
             obj.loadCompIM_powerPercent = loadCompIM_powerPercent;
             % Calculate IM params
             % (R)
             obj.loadCompIM_R = zeros(4,4,obj.numLoadComp);
             obj.loadCompIM_R(1,1,:) = loadCompIM_Rs;
             obj.loadCompIM_R(2,2,:) = loadCompIM_Rs;
             obj.loadCompIM_R(3,3,:) = loadCompIM_Rr;
             obj.loadCompIM_R(4,4,:) = loadCompIM_Rr;
             % (L)
             obj.loadCompIM_L = zeros(4,4,obj.numLoadComp);
             obj.loadCompIM_L(1,1,:) = loadCompIM_Lss;
             obj.loadCompIM_L(2,2,:) = loadCompIM_Lss;
             obj.loadCompIM_L(3,3,:) = loadCompIM_Lrr;
             obj.loadCompIM_L(4,4,:) = loadCompIM_Lrr;
             obj.loadCompIM_L(1,3,:) = loadCompIM_Lm;
             obj.loadCompIM_L(2,4,:) = loadCompIM_Lm;
             obj.loadCompIM_L(3,1,:) = loadCompIM_Lm;
             obj.loadCompIM_L(4,2,:) = loadCompIM_Lm;
             % (Swing)
             matSelect = [0 0; 0 0; 0 -1; 1 0]*[0 0 1 0;0 0 0 1];
             obj.loadCompIM_LTe = pagemtimes( ...
                 pagetranspose(obj.loadCompIM_L),matSelect);
             obj.loadCompIM_m = loadCompIM_m;
             obj.loadCompIM_H = loadCompIM_H;
             % Calculate ZIP
             assert(isequal(size(loadCompZip_PQabc),[2,obj.numLoadComp,3]))
             assert(all(sum(loadCompZip_PQabc,3)==1,'all'))
             obj.loadCompZip_PQabc = loadCompZip_PQabc;
             % ZIP filter for Vm
             assert(all(loadCompZip_Tau>0))
             obj.loadCompZip_Tau = loadCompZip_Tau;

         end

        % ***** GENERATOR PARAMS *****
        function obj = setGenrouParams(obj,genrou_bus, genrou_EsBase,genrou_IsBase, ...
                genrou_BaseMVA,...
                genrou_MW,genrou_Mvar, ...
                genrou_Ra,genrou_Rfd,genrou_R1d,genrou_R1q,genrou_R2q, ...
                genrou_Lad,genrou_Laq,genrou_L0,genrou_Ll, ...
                genrou_Lffd,genrou_Lf1d,genrou_L11d,genrou_L11q,genrou_L22q, ...
                genrou_H,genrou_KD)
            % Set bus numbers
            obj.genrou_bus = genrou_bus;
            obj.numGenrou = numel(genrou_bus);
            % Set base
            obj.genrou_EsBase = genrou_EsBase;
            obj.genrou_IsBase = genrou_IsBase;
            obj.genrou_IsBaseRMS = genrou_IsBase/sqrt(2);
            obj.genrou_BaseMVA = genrou_BaseMVA;
            assert(all(obj.genrou_BaseMVA~=0))
            % Set power
            obj.genrou_MW = genrou_MW;
            obj.genrou_Mvar = genrou_Mvar;
            % Set parameters
            obj.genrou_H = genrou_H;
            obj.genrou_KD = genrou_KD;
            % Winding resistance / inductance matrices, one 7x7 page
            % per machine
            [obj.genrou_RMatrix,obj.genrou_LMatrix,obj.genrou_LspMatrix] = ...
                setGenrouParams_matrices( ...
                genrou_Ra,genrou_Rfd,genrou_R1d,genrou_R1q,genrou_R2q, ...
                genrou_Lad,genrou_Laq,genrou_Ll,genrou_L0, ...
                genrou_Lffd,genrou_Lf1d,genrou_L11d,genrou_L11q,genrou_L22q);

            % For exc
            obj.genrou_Lad = genrou_Lad;
            obj.genrou_Rfd = genrou_Rfd;

            %
            obj.genrou_Xq = -obj.genrou_LMatrix(3:end,4);
            obj.genrou_Ra = obj.genrou_RMatrix(3:end,1);
            obj.genrou_Xd = -obj.genrou_LMatrix(3:end,1);
        end

        % ***** EXCITER PARAMS *****
        function obj = setExcDc1Params(obj,excDc1_ExcBus, excDc1_Exdc1Tr, ...
                excDc1_Exdc1Ta,excDc1_Exdc1Tc,excDc1_Exdc1Tb,excDc1_Exdc1Te, ...
                excDc1_Exdc1Tf1,excDc1_Exdc1Kf1,excDc1_Exdc1Ka,excDc1_Exdc1Ke)
            % Set parameters
            obj.excDc1_bus = excDc1_ExcBus;
            obj.numExcDc1 = numel(excDc1_ExcBus);

            obj.excDc1_Tr = excDc1_Exdc1Tr;
            obj.excDc1_Ta = excDc1_Exdc1Ta;
            obj.excDc1_Tc = excDc1_Exdc1Tc;
            obj.excDc1_Tb = excDc1_Exdc1Tb;
            obj.excDc1_Te = excDc1_Exdc1Te;
            obj.excDc1_Tf1 = excDc1_Exdc1Tf1;
            obj.excDc1_Kf1 = excDc1_Exdc1Kf1;
            obj.excDc1_Ka = excDc1_Exdc1Ka;
            obj.excDc1_Ke = excDc1_Exdc1Ke;

            [~,idxGenIsExc] = ismembertol(obj.excDc1_bus,obj.genrou_bus);
            obj.excDc1_genToExc = obj.genrou_Lad(idxGenIsExc)./obj.genrou_Rfd(idxGenIsExc);

        end
        % ***** GOV PARAMS *****
        function obj = setTgov1dParams(obj,tgov1d_Tgov1Bus, tgov1d_Tgov1T1,tgov1d_Tgov1T2, ...
                tgov1d_Tgov1T3,tgov1d_Tgov1R,tgov1d_Tgov1Dt)
            % Set parameters
            obj.tgov1d_bus = tgov1d_Tgov1Bus;
            obj.numTgov1d = numel(tgov1d_Tgov1Bus);

            obj.tgov1d_T1 = tgov1d_Tgov1T1;
            obj.tgov1d_T2 = tgov1d_Tgov1T2;
            obj.tgov1d_T3 = tgov1d_Tgov1T3;
            obj.tgov1d_R = tgov1d_Tgov1R; % .* (obj.baseMVA ./ obj.genrou_BaseMVA);
            obj.tgov1d_Dt = tgov1d_Tgov1Dt;

        end

        % *** IBR PARAMS ***
        function obj = setIbrGFMParams(obj,ibr_bus, ...
                ... % Base
                ibrGFM_BaseKv,ibrGFM_IBase,ibrGFM_MVA_base, ...
                ... % Power
                ibrGFM_MW,ibrGFM_MVar,...
                ... % Filter
                ibrGFM_Lf,ibrGFM_Rf,ibrGFM_Cf,ibrGFM_Rcap, ...
                ... % Coupling
                ibrGFM_Lc,ibrGFM_Rc, ...
                ... % Current controller
                ibrGFM_kC_i,ibrGFM_kC_p,ibrGFM_GC, ...
                ... % Voltage controller
                ibrGFM_kV_i,ibrGFM_kV_p,ibrGFM_GV, ...
                ... % Power
                ibrGFM_wmeas, ...
                ... % Droop
                ibrGFM_droopPercentP,ibrGFM_droopPercentQ, ...
                ... % Virtual resistor
                ibrGFM_w_VR,ibrGFM_Vr)

            % Bus
            obj.ibrGFM_bus = ibr_bus;
            obj.numIbrGFM = numel(obj.ibrGFM_bus);

            % Base
            obj.ibrGFM_BaseKV = ibrGFM_BaseKv;
            obj.ibrGFM_BaseKA = ibrGFM_IBase;
            obj.ibrGFM_MVA_base = ibrGFM_MVA_base;
            % Check there is factor of sqrt(3) because assumes VLL,RMS
            assert(all(abs(ibrGFM_BaseKv.*ibrGFM_IBase*sqrt(3)-ibrGFM_MVA_base)<obj.absTolCheck))

            % Power
            obj.ibrGFM_MW = ibrGFM_MW;
            obj.ibrGFM_MVar = ibrGFM_MVar;

            % Filter
            obj.ibrGFM_Lf = ibrGFM_Lf;
            obj.ibrGFM_Rf = ibrGFM_Rf;
            obj.ibrGFM_Cf = ibrGFM_Cf;
            obj.ibrGFM_Rcap = ibrGFM_Rcap;

            % Coupling
            obj.ibrGFM_Lc = ibrGFM_Lc;
            obj.ibrGFM_Rc = ibrGFM_Rc;

            % Current controller
            obj.ibrGFM_kC_i = ibrGFM_kC_i;
            obj.ibrGFM_kC_p = ibrGFM_kC_p;
            obj.ibrGFM_GC = ibrGFM_GC;
            % Voltage controller
            obj.ibrGFM_kV_i = ibrGFM_kV_i;
            obj.ibrGFM_kV_p = ibrGFM_kV_p;
            obj.ibrGFM_GV = ibrGFM_GV;
            % Power (LPF)
            obj.ibrGFM_wmeas = ibrGFM_wmeas;
            obj.ibrGFM_droopPercentP = ibrGFM_droopPercentP;
            obj.ibrGFM_droopPercentQ = ibrGFM_droopPercentQ;
            % Droop
            obj.ibrGFM_Mp = obj.ibrGFM_droopPercentP*obj.w0;
            obj.ibrGFM_Mq = obj.ibrGFM_droopPercentQ;
            % Virtual resistor
            obj.ibrGFM_w_VR = ibrGFM_w_VR;
            obj.ibrGFM_Rv = ibrGFM_Vr;

        end


        % *** PV PARAMS ***
        % *** LOAD ZIP + deltaF PARAMS ***
        % ***** node at C PARAMS *****

         %% Helper function for indexing

        % Gmatrix structure:
        % | [ bus - phA ] |
        % | [ bus - phB ] |
        % | [ bus - phC ] |
        % |---------------|
        % |---------------|
        % | [   gen1    ] | } gen # eqns
        % |      ...      |
        % | [   genx    ] |
        % |---------------|
        % | [   exc1    ] | } exc # eqns
        % | [    ...    ] |
        % | [   excx    ] |
        % |---------------|
        % | [   gov1    ] | } gove # eqns
        % | [    ...    ] |
        % | [   govx    ] |
        % |---------------|
        % | [   gov1    ] | } gov # eqns
        % | [    ...    ] |
        % | [   govx    ] |
        % |---------------|
        % |---------------|
        % | [  ibrGF1   ] | } grid forming IBR # eqns
        % | [    ...    ] |
        % | [  ibrGFx   ] |
        % |---------------|
        % |---------------|
        % | [ loadComp  ] | } composite load # eqns
        % | [    ...    ] |
        % | [ loadComp  ] |
        % |---------------|
        % | [  uel; oel ] | } exciter Q-limiting
        % | [    ...    ] |
        % | [  uel; oel  ] |


        % System
        function rowIdxAbc = getRowIdxAbc_fromBusNums(obj,busNum)
            % Must use "ismembertol" because multiple bus numbers
            [~,bus_idx] = ismembertol(busNum,obj.bus);
            rowIdxAbc = [bus_idx; obj.numBus+bus_idx; 2*obj.numBus+bus_idx];
        end
        function rowIdxAbc = getRowIdxAbc_fromBusNum(obj,busNum)
            % Use "==" and "find" instead of "ismembertol" because faster
            bus_idx = find(busNum==obj.bus);
            rowIdxAbc = bus_idx + [0;obj.numBus;2*obj.numBus];
        end
        function stampIdx = getStampIdx_static(obj)
            stampIdx = 1:obj.static_numStampVals;
        end
        % Gen
        function rowIdx = getRowIdx_fromGenrouIndices(obj,genrou_idx)
            rowIdx = 3*obj.numBus + 1+obj.numGenrouEqns*(genrou_idx-1);
        end
        function stampIdx = getStampIdx_fromGenrouIndices(obj,genrou_idx)
            stampIdx = obj.static_numStampVals + ...
                1+obj.genrou_numStampVals*(genrou_idx-1);
        end
        % Exc
        function rowIdx = getRowIdx_fromExcDc1Indices(obj,excDc1_idx)
            rowIdx = 3*obj.numBus + obj.numGenrouEqns*obj.numGenrou + ...
                1+obj.numExcDc1Eqns*(excDc1_idx-1);
        end
        function stampIdx = getStampIdx_fromExcDc1Indices(obj,excDc1_idx)
            stampIdx = obj.static_numStampVals + ...
                obj.genrou_numStampVals*obj.numGenrou + ...
                1+obj.excDc1_numStampVals*(excDc1_idx-1);
        end
        % Gov
        function rowIdx = getRowIdx_fromTgov1dIndices(obj,tgov1d_idx)
            rowIdx = 3*obj.numBus + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + 1+obj.numTgov1dEqns*(tgov1d_idx-1);
        end
        function stampIdx = getStampIdx_fromTgov1dIndices(obj,tgov1d_idx)
            stampIdx = obj.static_numStampVals + ...
                obj.genrou_numStampVals*obj.numGenrou + ...
                obj.excDc1_numStampVals*obj.numExcDc1 + ...
                1+obj.tgov1d_numStampVals*(tgov1d_idx-1);
        end
        % Grid forming IBR
        function rowIdx = getRowIdx_fromIbrGFMIndices(obj,ibrGFM_idx)
            rowIdx = 3*obj.numBus + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + 1+obj.numTgov1dEqns*obj.numTgov1d + ...
                obj.numIbrGFMEqns*(ibrGFM_idx-1);
        end
        function stampIdx = getStampIdx_fromIbrGFMIndices(obj,ibrGFM_idx)
            stampIdx = obj.static_numStampVals + ...
                obj.genrou_numStampVals*obj.numGenrou + ...
                obj.excDc1_numStampVals*obj.numExcDc1 + ...
                obj.tgov1d_numStampVals*obj.numTgov1d + ...
                1+obj.ibrGFM_numStampVals*(ibrGFM_idx-1);
        end
        % Composite loads
        function rowIdx = getRowIdx_fromLoadCompIndices(obj,loadComp_idx)
            rowIdx = 3*obj.numBus + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + 1+obj.numTgov1dEqns*obj.numTgov1d + ...
                obj.numIbrGFMEqns*obj.numIbrGFM + ...
                obj.numLoadCompEqns*(loadComp_idx-1);
        end
        function stampIdx = getStampIdx_fromLoadCompIndices(obj,loadComp_idx)
            stampIdx = obj.static_numStampVals + ...
                obj.genrou_numStampVals*obj.numGenrou + ...
                obj.excDc1_numStampVals*obj.numExcDc1 + ...
                obj.tgov1d_numStampVals*obj.numTgov1d + ...
                obj.ibrGFM_numStampVals*obj.numIbrGFM + ...
                1+obj.loadComp_numStampVals*(loadComp_idx-1);
        end
        % PV
        % Load ZIP
        % UEL
        % UEL
        % System frequency
        function rowIdx = getRowIdx_wSys(obj)
            rowIdx = obj.getSysSize;
        end
        % System
        function sysSize = getSysSize(obj)
            sysSize = 3*obj.numBus + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + obj.numTgov1dEqns*obj.numTgov1d + ...
                obj.numIbrGFMEqns*obj.numIbrGFM + ...
                obj.numLoadCompEqns*obj.numLoadComp + ...
                1;
        end
        function stampSize = getStampSize(obj)
            stampSize = obj.static_numStampVals + ...
                obj.genrou_numStampVals*obj.numGenrou + ...
                obj.excDc1_numStampVals*obj.numExcDc1 + ...
                obj.tgov1d_numStampVals*obj.numTgov1d + ...
                obj.ibrGFM_numStampVals*obj.numIbrGFM + ...
                obj.loadComp_numStampVals*obj.numLoadComp;
            % NOTE: The system-frequency row is part of the static stamps
            % already, so it needs no separate term here.
        end
        % function numStampVals = getNumSysFreqStampVals(obj)
        %     numStampVals = obj.numGenrou + obj.numIbrGFM + 1;
        % end
        % ^^ Don't need system frequency stamp explicitly since it is
        % calculated during linear stamping

        %% Helper function for PU (per-unit) conversions
        % Voltage base conventions:
        % - Bus   line-to-line, RMS, kV
        % - Gen   line-to-neutral, peak, V
        % - IBR   line-to-line, RMS, kV and kA
        % - Composite load IM   line-to-neutral, peak, V (same derivation
        %                       as the synchronous machine)

        % For current, use genrou_IsBaseRMS because need to convert between
        % RMS and peak when converting between bus and gen (NOTE: Do not
        % need to worry about line-to-line vs line-to-neutral because only
        % one type of current)
        function c = busToGenrou_I(obj,bus_idx,genrou_idx)
            c = obj.bus_IBase(bus_idx)*1000./obj.genrou_IsBaseRMS(genrou_idx);
        end
        function c = genrouToBus_I(obj,bus_idx,genrou_idx)
            c = obj.genrou_IsBaseRMS(genrou_idx)./(obj.bus_IBase(bus_idx)*1000);
        end
        %
        function c = busToIbrGFM_I(obj,bus_idx,ibrGFM_idx)
            c = obj.bus_IBase(bus_idx)./obj.ibrGFM_BaseKA(ibrGFM_idx);
        end
        function c = ibrGFMToBus_I(obj,bus_idx,ibrGFM_idx)
            c = obj.ibrGFM_BaseKA(ibrGFM_idx)./obj.bus_IBase(bus_idx);
        end
        %
        function c = busToLoadComp_I(obj,bus_idx,loadComp_idx)
            c = obj.bus_IBase(bus_idx)*1000./obj.loadComp_IsBaseRMS( ...
                loadComp_idx);
        end
        function c = loadCompToBus_I(obj,bus_idx,loadComp_idx)
            c = obj.loadComp_IsBaseRMS(loadComp_idx)./ ...
                (obj.bus_IBase(bus_idx)*1000);
        end

        % For voltage, the calculation is:
        % Vpu_bus = sqrt(3/2)*Vpu_gen*Vbase_gen/Vbase_bus
        % [LL,RMS]  [    LL,RMS     ]
        % ...but because our procedure uses:
        % Vbase_gen = Vbase_bus*sqrt(2/3)
        % ...we get:
        % Vpu_bus = sqrt(3/2)*Vpu_gen*( Vbase_bus*sqrt(2/3) )/Vbase_bus
        %         = sqrt(3/2)*sqrt(2/3)*Vpu_gen*Vbase_bus/Vbase_bus
        %         = Vpu_gen
        % ...so the scaling factor is just 1 when converting between the
        % two
        function c = busToGenrou_V(~,~,~)
             c = 1;
        end
        function c = genrouToBus_V(~,~,~)
            c = 1;
        end
        % For IBR, need conversion because not defining Vpu_IBR w.r.t bus
        % Vpu_bus = Vpu_ibr*Vbase_ibr/Vbase_bus
        % [LL,RMS]  [LL,RMS]
        function c = busToIbrGFM_V(obj,bus_idx,ibrGFM_idx)
             c = obj.bus_BaseKv(bus_idx)./obj.ibrGFM_BaseKV(ibrGFM_idx);
        end
        function c = ibrGFMToBus_V(obj,bus_idx,ibrGFM_idx)
            c = obj.ibrGFM_BaseKV(ibrGFM_idx)./obj.bus_BaseKv(bus_idx);
        end
        %
        function c = busToLoadComp_V(~,~,~)
            c = 1;
        end
        function c = loadCompToBus_V(~,~,~)
            c = 1;
        end


        %% Timer helper functions for OEL
    end

    methods(Static)
        %% Get device stamps
        % ***** TXLINE STAMPS *****
        function GTx = stampOfTxGMat(ReqBr,ReqSh)
            % Assume 3ph PI-line model (no coupling) as equiv Norton
            % Eqn @ "from" node:
            % Ifout + [(Vf-Vt)/ReqBr+IhistBr] + [Vf/ReqSh+IhistShF] = 0
            % Eqn @ "to" node:
            % Itout + [(Vt-Vf)/ReqBr-IhistBr] + [Vt/ReqSh+IhistShT] = 0
            % NOTE: x-vector is [Vf_a Vt_a|Vf_b Vt_b|Vf_c Vt_c]^T

            % Build conductance matrix
            GTx1ph = sparse([1 1 2 2],[1 2 1 2], ...
                [1/ReqBr+1/ReqSh -1/ReqBr -1/ReqBr 1/ReqBr+1/ReqSh]);
            % Convert to 3ph
            GTx = blkdiag(GTx1ph,GTx1ph,GTx1ph);
        end

        function [IHistBr,IHistShFrom,IHistShTo] = stampOfTxIhist( ...
                IBr_prev,IShFrom_prev,IShTo_prev, ...
                VFrom_prev,VTo_prev, ...
                aHistBr,bHistBr,aHistSh,bHistSh)
            % Branch
            IHistBr = EMT.stampOfIhist(IBr_prev,VFrom_prev,VTo_prev, ...
                aHistBr,bHistBr);
            % Shunt at "from" node (Assume currrent goes from node to gnd)
            IHistShFrom = EMT.stampOfIhist(IShFrom_prev,VFrom_prev,0, ...
                aHistSh,bHistSh);
            % Shunt at "to" node (Assume currrent goes from node to gnd)
            IHistShTo = EMT.stampOfIhist(IShTo_prev,VTo_prev,0, ...
                aHistSh,bHistSh);
        end

        % ***** XFMR STAMPS *****
        function GXfmr = stampOfXfmrGMat(Req)
            GXfmr = EMT.stampOfGMat(Req);
        end
        function IHistXfmr = stampOfXfmrIhist(I_prev,VFrom_prev,VTo_prev, ...
                aHist,bHist)
            IHistXfmr = EMT.stampOfIhist(I_prev,VFrom_prev,VTo_prev, ...
                aHist,bHist);
        end

        % ***** CONST LOAD STAMPS *****
        function GLoad = stampOfLoadGMat(Req)
            GLoad_full = EMT.stampOfGMat(Req);
            GLoad = GLoad_full([1 3 5],[1 3 5]);
        end
        function IHistLoad = stampOfLoadIhist(I_prev,VFrom_prev, ...
                aHist,bHist)
            IHistLoad = EMT.stampOfIhist(I_prev,VFrom_prev,0, ...
                aHist,bHist);
        end

        % ***** NODE C STAMPS *****
    end

    methods (Static)
        %% Helper functions that equivalent circuit stamps
        function GMat = stampOfGMat(Req)
            % Assume 3ph RL model (no coupling), and solve as equiv Norton
            % Eqn @ "from" node:
            % Ifout + (Vf-Vt)/Req + Ihist = 0
            % Eqn @ "to" node:
            % Itout + (Vt-Vf)/Req - Ihist = 0
            % NOTE: x-vector is [Vf_a Vt_a|Vf_b Vt_b|Vf_c Vt_c]^T

            % Build conductance matrix
            GMat1ph = sparse([1 1 2 2],[1 2 1 2], ...
                [1/Req -1/Req -1/Req 1/Req]);
            % Convert to 3ph
            GMat = blkdiag(GMat1ph,GMat1ph,GMat1ph);
        end
        function IHist = stampOfIhist(I_prev,VFrom_prev,VTo_prev, ...
                aHist,bHist)
            % NOTE: All IHist's and V's (inputs and outputs) are 3ph
            % NOTE: Assume current source points from "from" node to
            %       "to" node
            IHist = aHist.*I_prev + bHist.*(VFrom_prev - VTo_prev);
        end
        function [Req,a,b] = calcReqHistParams_RLSeriesWithDampingRp(...
                R,L,Rp,deltaT)
            % R - (L // Rp)
            % Equivalent resistance and history multipliers, with
            % Ihist = a*I(t-deltaT) + b*V(t-deltaT)
            invRp = 1./Rp;
            Req = (1 + R.*(deltaT./(2*L) + invRp)) ./ ...
                (deltaT./(2*L) + invRp);
            a = (1 - R.*(deltaT./(2*L) - invRp)) ./ ...
                (1 + R.*(deltaT./(2*L) + invRp));
            b = (deltaT./(2*L) - invRp) ./ ...
                (1 + R.*(deltaT./(2*L) + invRp));
        end
        function [Req,a,b, Req_C,a_C,b_C,c] =  ...
                calcReqHistParams_RLCSeriesWithDampingRp(...
                R,L,Rp,C,Rs,deltaT)
            % Branch topology:
            %   R - (L // Rp) - (C - Rs)
            % In general,  Req in parallel with current src Ihist:
            %   I(t) = V(t)/Req + Ihist(t-deltaT)
            %   Ihist(t-deltaT) = a*I(t-deltaT) + b*V(t-deltaT) +
            %                     c*VCseries(t-deltaT)
            %
            % Derivation:
            % Given:
            % I1 = V1/Req1 + Ihist1
            % I2 = V2/Req2 + Ihist
            % s.t. I1 = I2
            %
            % FIRST, LET'S SOLVE FOR REQ:
            % I1 = (V1+V2)/Req + Ihist
            %     1. Re-write OG equations:
            %        V1 = (I1 - Ihist1)*Req1
            %        V2 = (I2 - Ihist2)*Req2
            %        =>
            %        V1 + V2 = (I1 - Ihist1)*Req1 + (I2 - Ihist2)*Req2
            %     2. Substitute I1 = I2:
            %        V1 + V2 = (I1 - Ihist1)*Req1 + (I1 - Ihist2)*Req2
            %                = I1*(Req1 + Req2) - (Ihist1*Req1 + Ihist2*Req2)
            %     3. Solve for I1:
            %        I1 = (V1+V2)/(Req1+Req2) + ...
            %             (Ihist1*Req1 + Ihist2*Req2)/(Req1+Req2)
            %        => Req = Req1+Req2
            % SECOND, LET'S SOLVE FOR a and b (lowercase means prev time):
            % Ihist = a*(v1+v2) + b*i1
            %     1. From before, we know:
            %        Ihist = (Ihist1*Req1 + Ihist2*Req2)/(Req1+Req2)
            %     2. Substitute in definition of Ihist's, I2=I1, and Req from
            %        previous step 3:
            %        a*i1 + b*(v1+v2) = [(a1*v1 + b1*(v1+v2))*Req1 +
            %                            (a2*v2 + b2*(v1+v2))*Req2]/Req
            %     3. Re-combine terms involving voltages and currents on RHS:
            %        [(a1*Req1)*v1+(a2*Req2)*v2 + (b1*Req1+b2*Req2)*(v1+v2)]/Req
            %        ==> b = (b1*Req1+b2*Req2)/Req
            %     4. That leaves:
            %        a*i1 = [(a1*Req1)*v1+(a2*Req2)*v2]/Req
            %     5. Let V2 be an extra state we store (for capacitor):
            %        a*i1 = (a1*Req1/Req)*i1 +
            %                    (a2*Req2/Req-a1*Req1/Req)*v2
            [Req_RL,a_RL,b_RL] = EMT.calcReqHistParams_RLSeriesWithDampingRp(...
                R,L,Rp,deltaT);
            [Req_C,a_C,b_C] = EMT.calcReqHistParams_CWithDampingRs(C,Rs,deltaT);
            Req = Req_RL + Req_C;
            a = (a_RL*Req_RL+a_C*Req_C)/Req;
            b = (b_RL*Req_RL/Req);
            c = (b_C*Req_C/Req-b_RL*Req_RL/Req);
        end
        function [Req,a,b] = calcReqHistParams_CWithDampingRs(C,Rs,deltaT)
            % C - Rs
            % Equivalent resistance and history multipliers, with
            % Ihist = a*V(t-deltaT) + b*I(t-deltaT)
            Req = Rs + deltaT./(2*C);
            a = (2*C.*Rs - deltaT) ./ ...
                (2*C.*Rs + deltaT);
            b = -2*C./(2*C.*Rs + deltaT);
        end
        function [Req,a,b] = calcReqHistParams_RParCWithDampingRs( ...
                R,C,Rs,deltaT)
            % R // (C - Rs)
            % Req_R//C = 1/(1/R+(2*C/deltaT) = 1/(1/R+1/(deltaT/(2*C)))
            %                                = 1/(1/R+1/(Req_C))
            % =>
            % Req_R//(C-Rs) = 1/(1/R+1/[Req_C-Rs]) =
            %                 1/(1/R+1/[Rs+deltaT/(2*C)])
            Req = 1./(1./R + 1./(Rs + deltaT./(2*C)));
            % a_R//C is -1, so plug in a_C-Rs
            a = - (deltaT./(2*C) - Rs)./(deltaT./(2*C) + Rs);
            COEFF = -(Rs-deltaT./(2*C)) ./ (Rs+deltaT./(2*C));
            b = COEFF./R - 1./(Rs+deltaT./(2*C));
        end
        %% Calculate damping resistor values
        % Damping resistors that keep the companion models of a pure
        % inductor or capacitor well conditioned.
        function Rp = getDampingR_ParL(L,deltaT)
            Rp = (40/3)*L/deltaT;
        end
        function Rs = getDampingR_SeriesC(C,deltaT)
            Rs = (3/40)*deltaT./C;
        end


    end
end