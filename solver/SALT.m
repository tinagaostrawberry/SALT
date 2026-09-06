classdef SALT < handle
%SALT  Steady-state (harmonic-balance) solver built on EMT device models.
%
%   salt = SALT(rawFile,deviceFile)
%   salt = SALT(rawFile,deviceFile,Name,Value,...)
%
%   rawFile     PSS/E .raw (or any MATPOWER-loadable .m/.mat) case used for
%               the power-flow solution and the network data.
%   deviceFile  Spreadsheet of dynamic-device parameters. One sheet per
%               model; column 1 holds parameter names, each further column
%               is one device. Recognised sheets and the model each must
%               declare in its "type" row:
%
%                   gen               GENROU
%                   exc               EXDC1
%                   gov               TGOV1
%                   gfm ibr generic   GFM_GENERIC
%                   zip load          ZIP
%                   im load kundur    IM_KUNDUR
%
%   Name-value options (all optional; outage options default to false):
%
%     maxDeltaT      Largest time step for the derived one-period sample
%                    grid.                                 (default 1e-5)
%     loadScale      Multiplies the composite-load impedance, so a smaller
%                    scale means a larger load and a lower frequency.
%                                                          (default 1)
%     lineOutage     false, or indices into the transmission-line list to
%                    remove.                               (default false)
%     genOutage      false, or bus numbers whose machine is removed.
%                                                          (default false)
%     loadOutage     false, or bus numbers whose load is removed.
%                                                          (default false)
%     newSlackBus    Bus number to promote to slack, used when the outage
%                    removes the original slack machine.   (default [])
%     rerunPfInit    Solve the power flow before initialising, so the
%                    base case is a converged operating point.
%                                                          (default true)
%
%   The unknown vector x_v is laid out as
%
%       [ bus + MNA voltages, phases a/b/c, each as (re,im) ]   <- linear
%       [ GENROU | EXDC1 | TGOV1D | GFM IBR | composite load ]  <- devices
%       [ system frequency w ]                                  <- 1 row
%
%   The Jacobian is never written element by element. Every stamp is
%   emitted as a (value,row,col) triplet - see BUILDSTAMPMAPS for the row
%   and column maps, which are built once, and STAMPLINEAR / STAMPNONLINEAR
%   for the values, which are recomputed per Newton iteration across all
%   devices of a kind at once.

    properties (Constant)
        % System values
        f0 = 60                       % Assume 60 Hz
        w0 = 2*pi*SALT.f0
        genrouP = 2                   % Assume simple 2 pole generator

        % Equations per device
        numGenrouEqns = 15
        numExcDc1Eqns = 6             % 5 states + 1 vt calculation
        numTgov1dEqns = 3             % 2 states + 1 pmech calculation
        numIbrGFMEqns = 33            % 8 Park transf + 6 passive +
                                      % 4 current control + 4 voltage control +
                                      % 5 droop + 6 virtual resistor
        numLoadCompEqns = 20          % IM: 4 stator (volt,curr) + 2 rotor
                                      % (curr) + 2 speeds (wr,wslip) +
                                      % 1 system theta + 6 iphABC (re,im)
                                      % ZIP: 2 PQ + 2 Idq + 1 Vm

        % Triplet values per device (see BUILDSTAMPMAPS)
        numGenrouStampVals   = 73     % 68 in-block + 4 park-V + 1 speed
        numExcDc1StampVals   = 17     % 10 in-block + 4 vt + 3 into gen efd
        numTgov1dStampVals   = 11     % 5 in-block + 3 wr + 3 into gen Tm
        numIbrGFMStampVals   = 183    % 178 in-block + 4 park-V + 1 speed
        numLoadCompStampVals = 107    % 99 in-block + 4 park-V + 4 w_sys

        % Real and imaginary part
        numReImEqns = 2

        % Sanity check tolerance
        absTolCheck = 1e-6
    end

    properties
        %% General
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
        bus_Zbase

        % Power-flow case the network was built from (kept for post-processing)
        mpc
        caseName                      % network file basename, used to name outputs
        pf_busVm
        pf_busVa

        % System (Y*x_v = I, so I is negative)
        Y
        I
        x_v
        x_v0

        % Slack bus (extra equation, since frequency is a variable)
        slack_imScalConst
        slack_busNum


        %% Triplet stamping
        % Row/col maps are built once by BUILDSTAMPMAPS; the value arrays
        % are overwritten every Newton iteration.
        stampRow                      % [nnz x 1] row index of each triplet
        stampCol                      % [nnz x 1] col index of each triplet
        stampVal                      % [nnz x 1] value of each triplet
        stampSlice                    % struct of index ranges into the above
        % Local (device-relative) row/col patterns, one per device model
        map                           % struct of per-model index patterns

        %% Static devices
        % Txline - PI line (Br = branch, Sh = shunt)
        tx_busFrom
        tx_busTo
        tx_R
        tx_L
        tx_C
        tx_Smax
        % Xfmr - R-L series
        xfmr_busFrom
        xfmr_busTo
        xfmr_mnaIdxMap = zeros(0,2)   % Maps xfmr device -> extra MNA var
        xfmr_numMnaVars = 0
        xfmr_R
        xfmr_L
        % Constant-impedance load
        load_bus
        load_R
        load_X
        load_Z

        %% GENROU
        genrou_bus
        genrou_EsBase
        genrou_IsBase
        genrou_IsBaseRMS
        genrou_BaseMVA
        genrou_ZsBase
        genrou_LsBase
        genrou_MW0
        genrou_Mvar0
        genrou_RMatrix
        genrou_LMatrix
        genrou_LspMatrix
        genrou_H
        genrou_KD
        genrou_TorqueM0
        genrou_efd0
        genrou_Lad
        genrou_Rfd
        % Raw winding parameters, retained so an EMT twin can be built
        % straight from a solved SALT object (see BUILDEMT)
        genrou_R1d
        genrou_R1q
        genrou_R2q
        genrou_Laq
        genrou_L0
        genrou_Ll
        genrou_Lffd
        genrou_Lf1d
        genrou_L11d
        genrou_L11q
        genrou_L22q
        genrou_Xd
        genrou_Xdp
        genrou_Xqp
        genrou_Xdpp
        genrou_Xqpp
        genrou_Xq
        genrou_Xl
        genrou_Tdp
        genrou_Tdpp
        genrou_Tqp
        genrou_Tq0pp
        genrou_Ra

        %% Controllers
        % Controller parameters are per unit on the MACHINE base, not the
        % system base, and are used without conversion.

        % EXDC1 exciter
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
        excDc1_vref
        excDc1_genrouToExcDc1
        excDc1_genrouIdx             % which GENROU each exciter drives

        % TGOV1D governor
        tgov1d_bus
        tgov1d_T1
        tgov1d_T2
        tgov1d_T3
        tgov1d_R
        tgov1d_Dt
        tgov1d_Pref
        tgov1d_genrouIdx             % which GENROU each governor drives

        %% Grid-forming IBR
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
        ibrGFM_voDroop_q_set

        %% Composite load (IM + ZIP)
        % One induction motor in parallel with one static ZIP load. The
        % motor share of the bus MW is set per load by powerPercent.
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

        %% Time-domain sampling
        tSim
        deltaT
        maxDeltaT
    end

    methods
        %% ============================ Class ==========================
        function obj = SALT(rawFile,deviceFile,opts)
            arguments
                rawFile     (1,1) string
                deviceFile  (1,1) string
                opts.maxDeltaT   (1,1) double  = 1e-5
                opts.loadScale   (1,1) double  = 1
                opts.lineOutage                = false
                opts.genOutage                 = false
                opts.loadOutage                = false
                opts.newSlackBus               = []
                opts.rerunPfInit     (1,1) logical = true
            end

            % Steady-state model-index globals, used throughout the
            % stamping code. The time-domain set is EMT's business and is
            % populated by its constructor.
            saltInitIndices;

            obj.maxDeltaT = opts.maxDeltaT;
            obj = obj.updateTimeSimSamples;

            % ---- Power flow / network ---------------------------------
            pf = SALT.readPowerFlowCase(rawFile,opts.rerunPfInit);
            [~,obj.caseName] = fileparts(rawFile);
            obj.baseMVA = pf.baseMVA;
            obj.bus     = pf.bus_num;
            obj.numBus  = numel(pf.bus_num);
            obj.mpc     = pf.mpc;
            obj.pf_busVm = pf.bus_Vm;
            obj.pf_busVa = pf.bus_Va;

            % Slack bus. The frequency variable needs one extra equation
            % that pins the slack angle, so store the imag/real ratio.
            [obj.slack_busNum,obj.slack_imScalConst] = SALT.getSlackParams(pf.mpc);

            % Bus bases (bus_basekV is RMS, line-to-line)
            obj.bus_BaseKv = pf.bus_basekV;
            obj.bus_IBase  = obj.baseMVA ./ (obj.bus_BaseKv*sqrt(3));   % kA
            obj.bus_Zbase  = (obj.bus_BaseKv.^2) ./ obj.baseMVA;

            % ---- Device data file -------------------------------------
            dev = SALT.readDeviceFile(deviceFile);

            % ---- Outage bookkeeping (all optional, default off) --------
            idxLineOut = SALT.outageToIndex(opts.lineOutage);
            busGenOut  = SALT.outageToIndex(opts.genOutage);
            busLoadOut = SALT.outageToIndex(opts.loadOutage);

            % ---- Transmission lines (3ph PI, no mutual coupling) -------
            tx_busFrom = pf.line_from;
            tx_busTo   = pf.line_to;
            tx_R       = real(pf.line_RX);
            X          = imag(pf.line_RX);
            assert(all(X>0),'Transmission-line reactance must be inductive')
            tx_L       = X/obj.w0;
            tx_C       = pf.line_chg/(2*obj.w0);
            txSmax     = pf.line_Smax;
            if ~isempty(idxLineOut)
                keep = true(size(tx_busFrom)); keep(idxLineOut) = false;
                tx_busFrom = tx_busFrom(keep); tx_busTo = tx_busTo(keep);
                tx_R = tx_R(keep); tx_L = tx_L(keep); tx_C = tx_C(keep);
                txSmax = txSmax(keep);
            end
            obj = obj.setTxParams(tx_busFrom,tx_busTo,tx_R,tx_L,tx_C);
            obj.tx_Smax = txSmax;

            % ---- Transformers (3ph R-L, no mutual coupling) ------------
            obj = obj.setXfmrParams(pf.xfmr_from,pf.xfmr_to, ...
                real(pf.xfmr_RX),imag(pf.xfmr_RX)/obj.w0);

            % ---- Composite loads --------------------------------------
            idxLoadSys = 1:numel(pf.load_bus);
            if ~isempty(busLoadOut)
                idxLoadSys(ismember(pf.load_bus,busLoadOut)) = [];
            end
            % A bus that draws reactive power but no real power is a
            % shunt reactor or capacitor, not a motor. The composite load
            % bases its machine rating on real power (IM_BaseMVA in
            % setLoadCompFromData), so a zero-P bus would give a zero MVA
            % base and the per-unit conversion would divide by it. Those
            % buses go through the constant-impedance path below instead,
            % which is what they physically are.
            idxLoadZ   = idxLoadSys(pf.load_MW(idxLoadSys) == 0);
            idxLoadSys = idxLoadSys(pf.load_MW(idxLoadSys) ~= 0);
            obj = obj.setLoadCompFromData(pf,dev,idxLoadSys,opts.loadScale);

            % ---- Constant impedances to ground -------------------------
            % Fixed bus shunts - capacitor banks and reactors - together
            % with any zero-real-power load bus. Cases with neither (the
            % 4-bus and 39-bus) leave this empty.
            constZ_bus = pf.shunt_bus(:);
            constZ_Z   = 1./pf.shunt_Y(:);
            if ~isempty(idxLoadZ)
                [~,iBusZ] = ismembertol(pf.load_bus(idxLoadZ),pf.bus_num);
                % Same convention as the composite load: loadScale divides
                % the power, so it multiplies the impedance.
                Sz = (pf.load_MW(idxLoadZ) - 1j*pf.load_Mvar(idxLoadZ)) ...
                     ./ (obj.baseMVA*opts.loadScale);
                constZ_bus = [constZ_bus; pf.load_bus(idxLoadZ(:))];
                constZ_Z   = [constZ_Z;   (pf.bus_Vm(iBusZ(:)).^2)./Sz(:)];
            end
            if ~isempty(constZ_bus)
                obj = obj.setLoadParams(constZ_bus,constZ_Z);
            end

            % ---- Machines, controllers and inverters -------------------
            % Every power-flow generator must be claimed by exactly one
            % dynamic model: a GENROU on the "gen" sheet or a grid-forming
            % inverter on the "gfm ibr generic" sheet.
            busIbr = dev.ibr.("Bus #");
            busSg  = dev.gen.("Bus #");
            assert(isempty(intersect(busSg,busIbr)), ...
                'A bus cannot host both a GENROU and a GFM IBR')
            assert(isequal(sort([busSg(:);busIbr(:)]),sort(pf.gen_bus(:))), ...
                ['Device file does not cover the power-flow generators. ' ...
                 'Expected buses: ' mat2str(sort(pf.gen_bus(:)).')])

            if ~isempty(busGenOut) && ~isempty(opts.newSlackBus)
                % Include new slack bus if generator outage was applied to
                % slack
                obj.slack_busNum = opts.newSlackBus;
                [~,iSlk] = ismember(obj.slack_busNum,pf.bus_num);
                vs = pf.bus_Vm(iSlk)*exp(1j*pf.bus_Va(iSlk));
                obj.slack_imScalConst = imag(vs)/real(vs);
            end

            obj = obj.setGenrouFromData(pf,dev,busGenOut);
            obj = obj.setExcDc1FromData(dev,busGenOut);
            obj = obj.setTgov1dFromData(dev,busGenOut);
            obj = obj.setIbrGFMFromData(pf,dev,busGenOut);

            % ---- Initial guess + triplet maps -------------------------
            obj = obj.systemValInitializations(pf.bus_Vm,pf.bus_Va);
            obj = obj.buildStampMaps;
        end
    end

    methods (Static, Hidden)
        %% ========================== File parsing ======================
        function pf = readPowerFlowCase(fileName,rerunPfInit)
            % Load a PSS/E raw (or MATPOWER) case and solve the power flow.
            [~,~,ext] = fileparts(fileName);
            if strcmpi(ext,'.raw')
                mpc = psse2mpc(char(fileName),'',0);
            else
                mpc = loadcase(char(fileName));
            end
            if rerunPfInit
                mpc = runpf(mpc,mpoption('verbose',0,'out.all',0));
            end
            define_constants;

            pf.baseMVA    = mpc.baseMVA;
            pf.mpc        = mpc;
            pf.bus_num    = mpc.bus(:,1);
            pf.bus_type   = mpc.bus(:,2);
            pf.bus_Vm     = mpc.bus(:,8);
            pf.bus_Va     = deg2rad(mpc.bus(:,9));
            pf.bus_basekV = mpc.bus(:,10);

            % Fixed bus shunts, in MW / MVAr absorbed at 1.0 pu on the
            % system base. A capacitor bank has Bs > 0. SALT stamps these
            % through its constant-impedance load path.
            idxSh        = (mpc.bus(:,5)~=0) | (mpc.bus(:,6)~=0);
            pf.shunt_bus = mpc.bus(idxSh,1);
            pf.shunt_Y   = complex(mpc.bus(idxSh,5),mpc.bus(idxSh,6)) ...
                           / pf.baseMVA;

            % Loads
            % A bus is a load if it draws real OR reactive power.
            % This used to require both to be nonzero, which silently
            % dropped any bus with P but no Q - and left that power
            % unserved, so the frequency settled away from nominal.
            % The PF project (scriptSetupPF.m) has always used OR.
            idxLd        = (mpc.bus(:,3)~=0) | (mpc.bus(:,4)~=0);
            pf.load_bus  = mpc.bus(idxLd,1);
            pf.load_MW   = mpc.bus(idxLd,3);
            pf.load_Mvar = mpc.bus(idxLd,4);

            % Generators
            pf.gen_bus      = mpc.gen(:,1);
            pf.gen_MW       = mpc.gen(:,2);
            pf.gen_Mvar     = mpc.gen(:,3);
            pf.gen_MVA_base = mpc.gen(:,7);
            assert(all(pf.gen_MW+1j*pf.gen_Mvar~=0))

            % Branches: ratio 0 marks a line, anything else a transformer
            idxTx        = (mpc.branch(:,9)==0);
            pf.line_from = mpc.branch(idxTx,1);
            pf.line_to   = mpc.branch(idxTx,2);
            pf.line_RX   = complex(mpc.branch(idxTx,3),mpc.branch(idxTx,4));
            pf.line_chg  = mpc.branch(idxTx,5);
            pf.line_Smax = mpc.branch(idxTx,RATE_A);
            assert(all(pf.line_Smax>0))

            idxXfmr = ~idxTx;
            if any(mpc.branch(idxXfmr,9)~=1)
                warning('SALT only supports transformer ratios of 1.')
            end
            pf.xfmr_from = mpc.branch(idxXfmr,1);
            pf.xfmr_to   = mpc.branch(idxXfmr,2);
            pf.xfmr_RX   = complex(mpc.branch(idxXfmr,3),mpc.branch(idxXfmr,4));
        end

        function [slack_busNum,slack_imScalConst] = getSlackParams(mpc)
            % Bus type 3 is the slack bus in the MATPOWER convention.
            iSlk = find(mpc.bus(:,2)==3);
            slack_busNum = mpc.bus(iSlk,1);
            Vslack = mpc.bus(iSlk,8)*exp(1j*deg2rad(mpc.bus(iSlk,9)));
            % Pin the angle by fixing the imag/real ratio, which is
            % invariant to magnitude: angle(V) == angle(V*M)
            assert(real(Vslack)~=0)
            slack_imScalConst = imag(Vslack)/real(Vslack);
        end

        function dev = readDeviceFile(fileName)
            % Read every recognised sheet and check the declared model type
            % matches the model SALT implements.
            expect = struct( ...
                'gen',      {{"gen",             "GENROU"}}, ...
                'exc',      {{"exc",             "EXDC1"}}, ...
                'gov',      {{"gov",             "TGOV1"}}, ...
                'ibr',      {{"gfm ibr generic", "GFM_GENERIC"}}, ...
                'zip',      {{"zip load",        "ZIP"}}, ...
                'im',       {{"im load kundur",  "IM_KUNDUR"}});

            present = sheetnames(fileName);
            fn = fieldnames(expect);
            for k = 1:numel(fn)
                sheet = expect.(fn{k}){1};
                model = expect.(fn{k}){2};
                if ~any(strcmpi(present,sheet))
                    dev.(fn{k}) = table();
                    continue
                end
                T = SALT.readSheetTable(fileName,sheet);
                if isempty(T)
                    dev.(fn{k}) = table();
                    continue
                end
                assert(ismember('type',T.Properties.VariableNames), ...
                    "Sheet '"+sheet+"' has no 'type' row")
                bad = ~ismember(string(T.type),[model,""]);
                assert(~any(bad), ...
                    "Sheet '"+sheet+"' declares type '"+ ...
                    strjoin(unique(string(T.type(bad))),"', '")+ ...
                    "' but SALT only implements '"+model+"'")
                dev.(fn{k}) = T(ismember(string(T.type),model),:);
            end
        end

        function T = readSheetTable(fileName,sheetName)
            % Sheets are stored transposed: column 1 is the parameter name,
            % every further column is one device.
            C = readcell(fileName,'Sheet',sheetName);
            if size(C,2) < 2, T = table(); return, end
            body = C(:,2:end);
            % Drop columns that are entirely empty
            keep = false(1,size(body,2));
            for c = 1:size(body,2)
                keep(c) = ~all(cellfun(@(v) any(ismissing(v)),body(:,c)));
            end
            body = body(:,keep);
            if isempty(body), T = table(); return, end
            T = cell2table(body.','VariableNames',C(:,1));
        end

        function idx = outageToIndex(spec)
            % Outage options accept false (the default, meaning none) or a
            % list of indices / bus numbers.
            if isempty(spec) || (islogical(spec) && isscalar(spec) && ~spec)
                idx = [];
            else
                idx = spec(:);
            end
        end
    end

    methods (Hidden)
        %% ================= Device data -> device params ================
        function obj = setGenrouFromData(obj,pf,dev,busGenOut)
            T = dev.gen;
            if isempty(T), return, end
            busAll = T.("Bus #");
            keep = ~ismember(busAll,busGenOut);
            T = T(keep,:);
            if isempty(T), return, end
            busGenrou = T.("Bus #");

            % Locate each machine in the power-flow generator list
            [ok,iPf] = ismember(busGenrou,pf.gen_bus);
            assert(all(ok),'GENROU bus missing from the power flow')

            % --- Bases -------------------------------------------------
            % bus_basekV is RMS and line-to-line; the machine base is peak
            % and line-to-neutral, hence the sqrt(2/3).
            [~,iBus] = ismembertol(busGenrou,pf.bus_num);
            EsBase = pf.bus_basekV(iBus)*1000*sqrt(2/3);        % pk, LN (V)
            % Ibase_RMS = Sbase/(Vbase_LL_RMS*sqrt(3)), and with a
            % line-to-neutral machine base Vbase_LL_RMS = sqrt(3)*Vbase_LN_RMS,
            % so Ibase_PK = Sbase/(Vbase_LN_RMS*3/2).
            MVAbase = pf.gen_MVA_base(iPf);
            IsBase  = MVAbase*1e6./(EsBase*3/2);                % (A)
            obj.genrou_ZsBase = EsBase./IsBase;
            obj.genrou_LsBase = obj.genrou_ZsBase*1000/obj.w0;  % (mH)

            % --- Circuit matrices from the GENROU data sheet -----------
            p = SALT.genrouMatFromSheet(T);

            % --- Push into the class ----------------------------------
            obj = obj.setGenrouParams(busGenrou,EsBase,IsBase,MVAbase, ...
                pf.gen_MW(iPf),pf.gen_Mvar(iPf), ...
                p.Ra,p.Rfd,p.R1d,p.R1q,p.R2q, ...
                p.Lad,p.Laq,p.L0,p.Ll, ...
                p.Lffd,p.Lf1d,p.L11d,p.L11q,p.L22q, ...
                p.H,p.KD, ...
                p.Xd,p.Xdp,p.Xqp,p.Xdpp,p.Xqpp,p.Xq,p.Xl, ...
                p.Td0p,p.Tq0p,p.Td0pp,p.Tq0pp);
        end

        function obj = setExcDc1FromData(obj,dev,busGenOut)
            T = dev.exc;
            if isempty(T), return, end
            T = T(~ismember(T.("Bus #"),busGenOut),:);
            % An exciter only exists if its machine does
            T = T(ismember(T.("Bus #"),obj.genrou_bus),:);
            if isempty(T), return, end
            obj = obj.setExcDc1Params(T.("Bus #"),T.Tr,T.Ta,T.Tc,T.Tb, ...
                T.Te,T.Tf1,T.Kf1,T.Ka,T.Ke);
        end

        function obj = setTgov1dFromData(obj,dev,busGenOut)
            T = dev.gov;
            if isempty(T), return, end
            T = T(~ismember(T.("Bus #"),busGenOut),:);
            T = T(ismember(T.("Bus #"),obj.genrou_bus),:);
            if isempty(T), return, end
            obj = obj.setTgov1dParams(T.("Bus #"),T.("T1(>0)(sec)"), ...
                T.("T2(sec)"),T.("T3(>0)(sec)"),T.R,T.Dt);
        end

        function obj = setIbrGFMFromData(obj,pf,dev,busGenOut)
            T = dev.ibr;
            if isempty(T), return, end
            T = T(~ismember(T.("Bus #"),busGenOut),:);
            if isempty(T), return, end
            busIbr = T.("Bus #");

            % The inverter inherits the dispatch of the power-flow machine
            % it replaces, so no per-unit conversion is needed here.
            [ok,iPf] = ismember(busIbr,pf.gen_bus);
            assert(all(ok),'GFM IBR bus missing from the power flow')

            kvBase  = T.("Base kV");
            mvaBase = T.("MVA base");
            % Base current follows from a line-to-line, RMS voltage base
            % (unlike the machine, which uses peak line-to-neutral).
            kaBase  = mvaBase./(kvBase*sqrt(3));

            obj = obj.setIbrGFMParams(busIbr,kvBase,kaBase,mvaBase, ...
                pf.gen_MW(iPf),pf.gen_Mvar(iPf), ...
                T.Lf,T.Rf,T.("Cf (uF)")*1e-6,T.Rcap,T.Lc,T.Rc, ...
                T.kC_i,T.kC_p,T.GC, ...
                T.kV_i,T.kV_p,T.GV, ...
                T.wmeas,T.droopPercentP,T.droopPercentQ,T.w_VR,T.Rv);
        end

        function obj = setLoadCompFromData(obj,pf,dev,idxLoadSys,loadScale)
            Tz = dev.zip;
            Ti = dev.im;
            % A case with load buses but no ZIP / IM sheet is a data error,
            % not a case without load: silently returning here would drop
            % every load from the network and still converge, at a
            % completely different operating point.
            if isempty(idxLoadSys)
                return
            end
            sheetNames = ["'zip load'","'im load kundur'"];
            missing = sheetNames([isempty(Tz),isempty(Ti)]);
            assert(isempty(missing), ...
                ['The device file has %d composite load bus(es) to model ' ...
                 'but no %s sheet.'], numel(idxLoadSys), ...
                 strjoin(missing,' or '))

            % The sheets must cover every surviving load bus; they may
            % also list buses removed by a load outage, which are dropped.
            busLd = pf.load_bus(idxLoadSys);
            [okZ,iz] = ismember(busLd,Tz.("Bus #"));
            [okI,ii] = ismember(busLd,Ti.("Bus #"));
            assert(all(okZ),['ZIP load sheet is missing bus(es): ' ...
                mat2str(busLd(~okZ).')])
            assert(all(okI),['IM load sheet is missing bus(es): ' ...
                mat2str(busLd(~okI).')])
            Tz = Tz(iz,:);
            Ti = Ti(ii,:);
            nLd = numel(busLd);

            % --- Induction motor ---------------------------------------
            % Reactances are already per-unit at w_PU = 1, so no division
            % by w is needed to get the inductances.
            Lm  = Ti.Xm;
            Lss = Ti.Xs + Ti.Xm;
            Lrr = Ti.Xr + Ti.Xm;

            % Motor MVA base, taken as the motor's apparent power at full
            % loading. The motor R and L are given in PU on this base, so
            % it has to be sized consistently. PF only sizes the base; the
            % motor's reactive power is set by the slip, not by PF.
            powerPercent = Ti.powerPercent;
            IM_BaseMVA   = pf.load_MW(idxLoadSys).*powerPercent./Ti.PF;

            % Same base derivation as the synchronous machine: peak,
            % line-to-neutral.
            [~,iBus] = ismembertol(busLd,pf.bus_num);
            VBase = pf.bus_basekV(iBus)*1000*sqrt(2/3);        % pk, LN (V)
            IBase = IM_BaseMVA*1e6./(VBase*3/2);               % (A)

            % --- ZIP ---------------------------------------------------
            % Third dimension is ordered (Z,I,P).
            PQabc = zeros(2,nLd,3);
            PQabc(1,:,1) = Tz.Pz; PQabc(1,:,2) = Tz.Pi; PQabc(1,:,3) = Tz.Pp;
            PQabc(2,:,1) = Tz.Qz; PQabc(2,:,2) = Tz.Qi; PQabc(2,:,3) = Tz.Qp;

            % --- Dispatch, with the optional load perturbation ----------
            % loadScale multiplies the impedance: a smaller scale means a
            % larger load and therefore a lower system frequency.
            MW   = pf.load_MW(idxLoadSys)/loadScale;
            MVar = pf.load_Mvar(idxLoadSys)/loadScale;

            obj = obj.setLoadCompParams(busLd, ...
                VBase,IBase,IM_BaseMVA, MW,MVar, ...
                powerPercent, Ti.Rs,Ti.Rr, Lss,Lrr,Lm, Ti.H, Ti.m, ...
                PQabc, Tz.Tau.');
        end
    end

    methods (Static, Hidden)
        function p = genrouMatFromSheet(T)
            %GENROUMATFROMSHEET  GENROU nameplate data -> circuit matrices.
            %
            % Inductances and resistances are normalised PU, except Ra
            % (corresponding to idq0), which is not normalised with respect
            % to time. The swing-equation quantities H and KD are not
            % normalised.
            p.Xd    = T.Xd;
            p.Xq    = T.Xq;
            p.Xdp   = T.("X'd");
            p.Xqp   = T.("X'q");
            p.Xdpp  = T.("X''d=X''q");
            p.Xqpp  = T.("X''d=X''q");
            p.Xl    = T.Xl;
            p.Td0p  = T.("T'd0");
            p.Td0pp = T.("T''d0");
            p.Tq0p  = T.("T'q0");
            p.Tq0pp = T.("T''q0");
            p.H     = T.H;
            p.KD    = T.D;
            p.Ra    = T.Ra;

            % INDUCTANCES - D-axis
            p.Lad  = p.Xd - p.Xl;
            p.Lffd = (p.Lad.^2) ./ (p.Xd - p.Xdp);
            p.Lf1d = p.Lad;
            p.L11d = (p.Lad.^2) ./ (p.Xd - p.Xdpp);
            % INDUCTANCES - Q-axis
            p.Laq  = p.Xq - p.Xl;
            p.L11q = (p.Laq.^2) ./ (p.Xq - p.Xqp);
            p.L22q = (p.Laq.^2) ./ (p.Xq - p.Xdpp);
            % INDUCTANCES - other
            p.Ll = p.Xl;
            p.L0 = p.Xl;                 % X0 is taken equal to Xl
            % RESISTANCES - D-axis
            Lfd = (p.Xdp - p.Xl).*p.Lad ./ (p.Lad - (p.Xdp - p.Xl));
            z   = p.Xdpp - p.Xl;
            y   = p.Lad.*Lfd ./ (p.Lad + Lfd);
            L1d = y.*z ./ (y - z);
            p.R1d = (y + L1d) ./ p.Td0pp;
            p.Rfd = (p.Lad + Lfd) ./ p.Td0p;
            % RESISTANCES - Q-axis
            L1q = (p.Xqp - p.Xl).*p.Laq ./ (p.Laq - (p.Xqp - p.Xl));
            z   = p.Xdpp - p.Xl;
            y   = p.Laq.*L1q ./ (p.Laq + L1q);
            L2q = y.*z ./ (y - z);
            p.R2q = (y + L2q) ./ p.Tq0pp;
            p.R1q = (p.Laq + L1q) ./ p.Tq0p;
        end
    end

    methods
        %% ======================= Set device params ====================

        % ***** TXLINE *****
        function obj = setTxParams(obj,tx_busFrom,tx_busTo, tx_R,tx_L,tx_C)
            obj.tx_busFrom = tx_busFrom(:);
            obj.tx_busTo   = tx_busTo(:);
            obj.tx_R = tx_R(:);
            obj.tx_L = tx_L(:);
            obj.tx_C = tx_C(:);
            assert(all(obj.tx_C~=0),'PI-line shunt capacitance must be nonzero')
            assert(all(obj.tx_R~=0), ...
                ['A lossless line needs an MNA branch current, which is ' ...
                 'not supported while frequency is a variable'])
        end

        % ***** XFMR *****
        function obj = setXfmrParams(obj,xfmr_busFrom,xfmr_busTo,xfmr_R,xfmr_L)
            obj.xfmr_busFrom = xfmr_busFrom(:);
            obj.xfmr_busTo   = xfmr_busTo(:);
            obj.xfmr_R = xfmr_R(:);
            obj.xfmr_L = xfmr_L(:);
            % A lossless transformer has no admittance, so it needs an
            % explicit branch-current unknown (modified nodal analysis).
            deviceIdx = find(obj.xfmr_R==0);
            if ~isempty(deviceIdx)
                n = numel(deviceIdx);
                varIdx = obj.xfmr_numMnaVars + (1:n).';
                obj.xfmr_mnaIdxMap = [obj.xfmr_mnaIdxMap; [deviceIdx varIdx]];
                obj.xfmr_numMnaVars = obj.xfmr_numMnaVars + n;
            end
        end

        % ***** CONSTANT-IMPEDANCE LOAD *****
        function obj = setLoadParams(obj,load_bus,load_Z)
            obj.load_bus = load_bus(:);
            obj.load_Z = load_Z(:);
            obj.load_X = imag(obj.load_Z);
            obj.load_R = real(obj.load_Z);
            % A shunt element to ground is stamped as a plain admittance:
            % even a lossless reactor has the finite admittance 1/(jwL),
            % and its derivative with respect to w is finite too, so no MNA
            % branch current is needed. (That is not true of a lossless
            % element in series between two buses, which is why the
            % transformer path does use MNA.)
            assert(all(isfinite(obj.load_Z)), ...
                'Shunt impedance must be finite')
            assert(~any((obj.load_R==0)&(obj.load_X==0)), ...
                'A shunt with zero impedance is a short circuit to ground')
        end

        % ***** GENROU *****
        function obj = setGenrouParams(obj,genrou_bus, ...
                genrou_EsBase,genrou_IsBase,genrou_BaseMVA, ...
                genrou_MW0,genrou_Mvar0, ...
                Ra,Rfd,R1d,R1q,R2q, ...
                Lad,Laq,L0,Ll, Lffd,Lf1d,L11d,L11q,L22q, ...
                H,KD, ...
                Xd,Xdp,Xqp,Xdpp,Xqpp,Xq,Xl, ...
                Td0p,Tq0p,Td0pp,Tq0pp)

            obj.genrou_bus = genrou_bus(:);
            obj.numGenrou  = numel(obj.genrou_bus);

            obj.genrou_EsBase    = genrou_EsBase(:);
            obj.genrou_IsBase    = genrou_IsBase(:);
            obj.genrou_IsBaseRMS = genrou_IsBase(:)/sqrt(2);
            obj.genrou_BaseMVA   = genrou_BaseMVA(:);

            % Not in machine PU nor bus PU
            obj.genrou_MW0   = genrou_MW0(:);
            obj.genrou_Mvar0 = genrou_Mvar0(:);

            obj.genrou_H  = H(:);
            obj.genrou_KD = KD(:);

            % Winding resistance / inductance matrices, one 7x7 page
            % per machine
            [obj.genrou_RMatrix,obj.genrou_LMatrix,obj.genrou_LspMatrix] = ...
                setGenrouParams_matrices( ...
                Ra,Rfd,R1d,R1q,R2q, Lad,Laq,Ll,L0, ...
                Lffd,Lf1d,L11d,L11q,L22q);

            % Needed by the exciter
            obj.genrou_Lad = Lad(:);
            obj.genrou_Rfd = Rfd(:);
            % Retained verbatim for BUILDEMT
            obj.genrou_R1d = R1d(:);   obj.genrou_R1q = R1q(:);
            obj.genrou_R2q = R2q(:);   obj.genrou_Laq = Laq(:);
            obj.genrou_L0  = L0(:);    obj.genrou_Ll  = Ll(:);
            obj.genrou_Lffd = Lffd(:); obj.genrou_Lf1d = Lf1d(:);
            obj.genrou_L11d = L11d(:); obj.genrou_L11q = L11q(:);
            obj.genrou_L22q = L22q(:);

            obj.genrou_Xd = Xd(:);   obj.genrou_Xdp = Xdp(:);
            obj.genrou_Xqp = Xqp(:); obj.genrou_Xdpp = Xdpp(:);
            obj.genrou_Xqpp = Xqpp(:); obj.genrou_Xq = Xq(:);
            obj.genrou_Xl = Xl(:);
            obj.genrou_Tdp = Td0p(:);  obj.genrou_Tdpp = Td0pp(:);
            obj.genrou_Tqp = Tq0p(:);  obj.genrou_Tq0pp = Tq0pp(:);
            obj.genrou_Ra = Ra(:);
        end

        % ***** EXDC1 *****
        function obj = setExcDc1Params(obj,bus,Tr,Ta,Tc,Tb,Te,Tf1,Kf1,Ka,Ke)
            obj.excDc1_bus = bus(:);
            obj.numExcDc1  = numel(obj.excDc1_bus);
            if obj.numExcDc1==0, return, end
            obj.excDc1_Tr = Tr(:);  obj.excDc1_Ta = Ta(:);
            obj.excDc1_Tc = Tc(:);  obj.excDc1_Tb = Tb(:);
            obj.excDc1_Te = Te(:);  obj.excDc1_Tf1 = Tf1(:);
            obj.excDc1_Kf1 = Kf1(:);
            obj.excDc1_Ka = Ka(:);  obj.excDc1_Ke = Ke(:);
            [~,iGen] = ismembertol(obj.excDc1_bus,obj.genrou_bus);
            obj.excDc1_genrouIdx = iGen;
            obj.excDc1_genrouToExcDc1 = ...
                obj.genrou_Lad(iGen)./obj.genrou_Rfd(iGen);
        end

        % ***** TGOV1D *****
        function obj = setTgov1dParams(obj,bus,T1,T2,T3,R,Dt)
            obj.tgov1d_bus = bus(:);
            obj.numTgov1d  = numel(obj.tgov1d_bus);
            if obj.numTgov1d==0, return, end
            obj.tgov1d_T1 = T1(:); obj.tgov1d_T2 = T2(:);
            obj.tgov1d_T3 = T3(:); obj.tgov1d_R = R(:);
            obj.tgov1d_Dt = Dt(:);
            [~,iGen] = ismembertol(obj.tgov1d_bus,obj.genrou_bus);
            obj.tgov1d_genrouIdx = iGen;
        end

        % ***** GRID-FORMING IBR *****
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
            % Filter
            obj.ibrGFM_Lf = Lf(:); obj.ibrGFM_Rf = Rf(:);
            obj.ibrGFM_Cf = Cf(:); obj.ibrGFM_Rcap = Rcap(:);
            % Coupling
            obj.ibrGFM_Lc = Lc(:); obj.ibrGFM_Rc = Rc(:);
            % Current controller
            obj.ibrGFM_kC_i = kC_i(:); obj.ibrGFM_kC_p = kC_p(:);
            obj.ibrGFM_GC = GC(:);
            % Voltage controller
            obj.ibrGFM_kV_i = kV_i(:); obj.ibrGFM_kV_p = kV_p(:);
            obj.ibrGFM_GV = GV(:);
            % Power measurement low-pass filter and droop
            obj.ibrGFM_wmeas = wmeas(:);
            obj.ibrGFM_droopPercentP = droopPercentP(:);
            obj.ibrGFM_droopPercentQ = droopPercentQ(:);
            % Virtual resistor
            obj.ibrGFM_w_VR = w_VR(:); obj.ibrGFM_Rv = Rv(:);
        end

        % ***** COMPOSITE LOAD *****
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
            % Retained verbatim for BUILDEMT
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
    end

    methods
        %% ==================== Indexing helpers ========================
        % ------------------------ x_v structure -----------------------
        %                   ___________      ____________  ---> BUS/LINEAR
        % |     BUS 1     | } Re,Im    |                 |
        % |     BUS 2     | } Re,Im    | PH A            |
        % |     ...       |            |                 |
        % |    MNA VARS   | } Re,Im    |                 |
        % |     ...       | ___________|                 |
        % |---------------| ___________                  |
        % |     BUS 1     | } Re,Im    |                 |
        % |     ...       |            | PH B            |
        % |    MNA VARS   | } Re,Im    |                 |
        % |     ...       | ___________|                 |
        % |---------------| ___________                  |
        % |     BUS 1     | } Re,Im    |                 |
        % |     ...       |            | PH C            |
        % |    MNA VARS   | } Re,Im    |                 |
        % |     ...       | ___________|      ___________|
        % |---------------| ___________       ___________ ---> NONLIN DEV
        % | NONLIN VAR 1  | } Re,Im    |                 |     In order:
        % | NONLIN VAR 2  | } Re,Im    | NONLIN #1       |     1. GENROU
        % |     ...       | ___________|                 |     2. EXDC1
        % |---------------| ___________                  |     3. TGOV1D
        % | NONLIN VAR 1  | } Re,Im    | NONLIN #2       |     4. GFM IBR
        % |     ...       | ___________|      ___________|     5. ZIP+IM
        % |---------------|                   ___________
        % |  SYSTEM FREQ  |                              | ---> 1 row
        % |---------------|                   ___________|

        function n = getNumMnaVars(obj)
            n = obj.xfmr_numMnaVars;
        end

        function rowIdxAbc = getRowIdxAbc_fromBusNums(obj,busNum)
            % Must use "ismembertol" because multiple bus numbers
            [~,bus_idx] = ismembertol(busNum,obj.bus);
            rowIdxAbc = obj.busIdxToRowAbc(bus_idx(:).');
        end
        function rowIdxAbc = getRowIdxAbc_fromBusNum(obj,busNum)
            % Use "==" and "find" instead of "ismembertol" because faster
            rowIdxAbc = obj.busIdxToRowAbc(find(busNum==obj.bus)); %#ok<*FNDSB>
        end
        function rowIdxAbc = busIdxToRowAbc(obj,bus_idx)
            % [N x 1] phase-A rows broadcast against the [1 x 3] phase
            % offsets, then flattened column-wise to a [3N x 1] vector.
            phaseA = (obj.numReImEqns*(bus_idx(:).'-1) + (1:obj.numReImEqns).');
            phaseOffset = (obj.numReImEqns*(obj.numBus+obj.getNumMnaVars)) ...
                * [0, 1, 2];
            rowIdxAbc = reshape(phaseA,[],1) + phaseOffset;
            rowIdxAbc = rowIdxAbc(:);
        end
        function rowIdxAbc = getRowIdxAbc_MNA_fromMnaVarIdx(obj,varIdx)
            rowOnPhase = obj.numReImEqns*obj.numBus + ...
                obj.numReImEqns*(reshape(varIdx,1,[])-1) + ...
                (1:obj.numReImEqns).';
            off = obj.numReImEqns*(obj.numBus+obj.getNumMnaVars)*[0,1,2];
            rowIdxAbc = reshape(rowOnPhase,[],1) + off;
            rowIdxAbc = rowIdxAbc(:);
        end

        % (Nonlinear device blocks)
        function rowIdx = getRowIdx_fromGenrouIndices(obj,idx)
            rowIdx = obj.getLinSysSize + 1 + obj.numGenrouEqns*(idx-1);
        end
        function rowIdx = getRowIdx_fromExcDc1Indices(obj,idx)
            rowIdx = obj.getLinSysSize + obj.numGenrouEqns*obj.numGenrou + ...
                1 + obj.numExcDc1Eqns*(idx-1);
        end
        function rowIdx = getRowIdx_fromTgov1dIndices(obj,idx)
            rowIdx = obj.getLinSysSize + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + ...
                1 + obj.numTgov1dEqns*(idx-1);
        end
        function rowIdx = getRowIdx_fromIbrGFMIndices(obj,idx)
            rowIdx = obj.getLinSysSize + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + ...
                obj.numTgov1dEqns*obj.numTgov1d + ...
                1 + obj.numIbrGFMEqns*(idx-1);
        end
        function rowIdx = getRowIdx_fromLoadCompIndices(obj,idx)
            rowIdx = obj.getLinSysSize + obj.numGenrouEqns*obj.numGenrou + ...
                obj.numExcDc1Eqns*obj.numExcDc1 + ...
                obj.numTgov1dEqns*obj.numTgov1d + ...
                obj.numIbrGFMEqns*obj.numIbrGFM + ...
                1 + obj.numLoadCompEqns*(idx-1);
        end
        % (Frequency variable - always the very last row)
        function rowIdx = getRowIdx_w(obj)
            rowIdx = obj.getSysSize;
            if ~isempty(obj.x_v)
                assert(rowIdx==length(obj.x_v))
            end
        end

        function linSysSize = getLinSysSize(obj)
            linSysSize = 3*obj.numReImEqns*(obj.numBus+obj.getNumMnaVars);
        end
        function sysSize = getSysSize(obj)
            sysSize = obj.getLinSysSize + ...
                obj.numGenrou*obj.numGenrouEqns + ...
                obj.numExcDc1*obj.numExcDc1Eqns + ...
                obj.numTgov1d*obj.numTgov1dEqns + ...
                obj.numIbrGFM*obj.numIbrGFMEqns + ...
                obj.numLoadComp*obj.numLoadCompEqns + ...
                1;                                  % system frequency
        end

        %% ================= Per-unit (PU) conversions ==================
        % Voltage base conventions:
        % - Bus                 line-to-line, RMS, kV
        % - GENROU              line-to-neutral, peak, V
        % - GFM IBR             line-to-line, RMS, kV and kA
        % - Composite-load IM   line-to-neutral, peak, V (same derivation
        %                       as the synchronous machine)

        % For current, use the RMS machine base because the bus is RMS and
        % the machine is peak. (Only one kind of current, so line-to-line
        % vs line-to-neutral does not enter.)
        function c = busToGenrou_I(obj,bus_idx,idx)
            c = obj.bus_IBase(bus_idx)*1000./obj.genrou_IsBaseRMS(idx);
        end
        function c = genrouToBus_I(obj,bus_idx,idx)
            c = obj.genrou_IsBaseRMS(idx)./(obj.bus_IBase(bus_idx)*1000);
        end
        function c = busToIbrGFM_I(obj,bus_idx,idx)
            c = obj.bus_IBase(bus_idx)./obj.ibrGFM_BaseKA(idx);
        end
        function c = ibrGFMToBus_I(obj,bus_idx,idx)
            c = obj.ibrGFM_BaseKA(idx)./obj.bus_IBase(bus_idx);
        end
        function c = busToLoadComp_I(obj,bus_idx,idx)
            c = obj.bus_IBase(bus_idx)*1000./obj.loadComp_IsBaseRMS(idx);
        end
        function c = loadCompToBus_I(obj,bus_idx,idx)
            c = obj.loadComp_IsBaseRMS(idx)./(obj.bus_IBase(bus_idx)*1000);
        end

        % For voltage the conversion is
        %   Vpu_bus = sqrt(3/2)*Vpu_gen*Vbase_gen/Vbase_bus
        %   [LL,RMS]           [    LL,RMS     ]
        % and since our procedure uses Vbase_gen = Vbase_bus*sqrt(2/3),
        %   Vpu_bus = sqrt(3/2)*sqrt(2/3)*Vpu_gen = Vpu_gen
        % so the machine and composite-load scale factors are just 1.
        function c = busToGenrou_V(~,~,~),    c = 1; end
        function c = genrouToBus_V(~,~,~),    c = 1; end
        function c = busToLoadComp_V(~,~,~),  c = 1; end
        function c = loadCompToBus_V(~,~,~),  c = 1; end
        % The IBR PU voltage is not referred to the bus base, so it does
        % need a conversion: Vpu_bus = Vpu_ibr*Vbase_ibr/Vbase_bus
        function c = busToIbrGFM_V(obj,bus_idx,idx)
            c = obj.bus_BaseKv(bus_idx)./obj.ibrGFM_BaseKV(idx);
        end
        function c = ibrGFMToBus_V(obj,bus_idx,idx)
            c = obj.ibrGFM_BaseKV(idx)./obj.bus_BaseKv(bus_idx);
        end

        %% ================= Time-domain sample grid ====================
        function obj = updateTimeSimSamples(obj)
            if isempty(obj.x_v)
                endTime = 1/obj.f0;
            else
                endTime = 1/(obj.x_v(obj.getRowIdx_w)/(2*pi));
            end
            obj.deltaT = endTime/ceil(endTime/obj.maxDeltaT);
            obj.tSim = (0:obj.deltaT:endTime).';
        end

        function x = getTimeDomainSignal_iDFT(obj,x_v_idx)
            % Linear rows oscillate at the system frequency; device states
            % are DC in the rotating frame.
            if x_v_idx <= obj.getLinSysSize
                freq = obj.x_v(obj.getRowIdx_w)/(2*pi);
            else
                freq = 0;
            end
            x = getTimeDomainSignal_iDFT(obj.tSim,@(~)x_v_idx,obj.x_v,freq,1);
        end
    end

    methods
        %% ==================== Initialization ==========================
        function obj = systemValInitializations(obj,bus_Vm,bus_Va)
            % Machine setpoints and internal states come from the
            % unperturbed operating point; load setpoints and internal
            % states are those of the perturbed system.
            global idxGenrou_Idq idxGenrou_Ifd idxGenrou_Edq idxGenrou_Efd ...
                   idxGenrou_Wr idxGenrou_Theta idxGenrou_IphABC idxGenrou_Tm %#ok<*GVMIS>

            obj.x_v = zeros(obj.getSysSize,1);

            % Bus voltages at the fundamental only.
            % NOTE: x_v bus voltage is line-to-line (LL), RMS
            [bus_VComplex_3ph,bus_VComplex] = busValInitializations(bus_Vm,bus_Va);
            idxBus = obj.getRowIdxAbc_fromBusNums(obj.bus);
            obj.x_v(idxBus(1:2:end)) = real(bus_VComplex_3ph);
            obj.x_v(idxBus(2:2:end)) = imag(bus_VComplex_3ph);

            % Machines (also fixes the initial efd and Tm setpoints)
            if obj.numGenrou > 0
                [obj,genrou_VAbsPhA,efd,Tm] = ...
                    obj.genrouValInitializations_assumeW0(bus_VComplex,bus_Va);
                obj.genrou_efd0 = efd;
                obj.genrou_TorqueM0 = Tm;
            end
            % Exciters
            if obj.numExcDc1 > 0
                et = genrou_VAbsPhA(obj.excDc1_genrouIdx);   % terminal voltage
                [obj,vref] = obj.excDc1ValInitializations_assumeW0(et);
                obj.excDc1_vref = vref;
            end
            % Governors
            if obj.numTgov1d > 0
                [obj,pref] = obj.tgov1dValInitializations_assumeW0;
                obj.tgov1d_Pref = pref;
            end
            % Grid-forming IBRs
            if obj.numIbrGFM > 0
                [obj,pset,qset,voDroop_q_set] = ...
                    obj.ibrGFMValInitializations_assumeW0(bus_VComplex,bus_Va);
                obj.ibrGFM_pset = pset;
                obj.ibrGFM_qset = qset;
                obj.ibrGFM_voDroop_q_set = voDroop_q_set;
                obj.ibrGFM_Mp = obj.ibrGFM_droopPercentP*obj.w0;
                obj.ibrGFM_Mq = obj.ibrGFM_droopPercentQ;
            end
            % Composite loads
            if obj.numLoadComp > 0
                [obj,Tm0,PQz,PQi,PQp] = ...
                    obj.loadCompValInitializations_assumeW0(bus_VComplex);
                obj.loadCompIM_Tm0 = Tm0;
                obj.loadCompZip_PQz = PQz;
                obj.loadCompZip_PQi = PQi;
                obj.loadCompZip_PQp = PQp;
            end

            % Transformer branch currents for the MNA unknowns
            for k = 1:obj.xfmr_numMnaVars
                deviceIdx = obj.xfmr_mnaIdxMap(k,1);
                mnaVarIdx = obj.xfmr_mnaIdxMap(k,2);
                idxF = obj.getRowIdxAbc_fromBusNum(obj.xfmr_busFrom(deviceIdx));
                idxT = obj.getRowIdxAbc_fromBusNum(obj.xfmr_busTo(deviceIdx));
                VF = complex(obj.x_v(idxF(1:2:end)),obj.x_v(idxF(2:2:end)));
                VT = complex(obj.x_v(idxT(1:2:end)),obj.x_v(idxT(2:2:end)));
                assert(obj.xfmr_R(deviceIdx)==0)
                Yeq = obj.calcYeq_RLSeries(0,obj.xfmr_L(deviceIdx),obj.w0);
                I_F_T = (VF-VT)*Yeq;
                idxI = obj.getRowIdxAbc_MNA_fromMnaVarIdx(mnaVarIdx);
                obj.x_v(idxI(1:2:end)) = real(I_F_T);
                obj.x_v(idxI(2:2:end)) = imag(I_F_T);
            end

            % Frequency
            obj.x_v(obj.getRowIdx_w) = obj.w0;
            obj.x_v0 = obj.x_v;

            % Silence the unused-global warnings for indices consumed by
            % the helper calls above
            assert(~isempty(idxGenrou_Idq) && ~isempty(idxGenrou_Ifd) && ...
                   ~isempty(idxGenrou_Edq) && ~isempty(idxGenrou_Efd) && ...
                   ~isempty(idxGenrou_Wr)  && ~isempty(idxGenrou_Theta) && ...
                   ~isempty(idxGenrou_IphABC) && ~isempty(idxGenrou_Tm))
        end

        %% ------------- Per-device initialization helpers --------------
        function [obj,genrou_VAbsPhA,genrou_efd,genrou_TorqueM] = ...
                genrouValInitializations_assumeW0(obj,bus_VComplexPhA,bus_Va)
            % Set i (1:7), e (8:12), wr and theta
            global idxGenrou_Idq idxGenrou_Ifd idxGenrou_Edq idxGenrou_Efd ...
                   idxGenrou_Wr idxGenrou_Theta idxGenrou_IphABC idxGenrou_Tm

            [id,iq,ifd,ed,eq,genrou_efd,THETA,genrou_TorqueM, ...
                genrou_VAbsPhA,VComplexPhA_busPU,IComplexPhA_busPU] = ...
                genrouValInitCalculations_assumeW0(obj,bus_VComplexPhA,bus_Va);

            % Only the real parts are set, since these states sit at DC
            r = obj.getRowIdx_fromGenrouIndices(1:obj.numGenrou)-1;
            obj.x_v(r+idxGenrou_Idq(1)) = id;
            obj.x_v(r+idxGenrou_Idq(2)) = iq;
            obj.x_v(r+idxGenrou_Ifd)    = ifd;
            obj.x_v(r+idxGenrou_Edq(1)) = ed;
            obj.x_v(r+idxGenrou_Edq(2)) = eq;
            obj.x_v(r+idxGenrou_Efd)    = genrou_efd;
            obj.x_v(r+idxGenrou_Wr)     = obj.w0;      % Assume nominal
            obj.x_v(r+idxGenrou_Theta)  = THETA;
            obj.x_v(r+idxGenrou_Tm)     = genrou_TorqueM;

            % Terminal current in abc, at nominal frequency
            for k = 1:obj.numGenrou
                [Ia,Ib,Ic] = calc3Ph(IComplexPhA_busPU(k));
                obj.x_v(r(k)+idxGenrou_IphABC) = ...
                    [real(Ia);imag(Ia);real(Ib);imag(Ib);real(Ic);imag(Ic)];
            end

            % Sanity check on the MW/MVar the initialization reproduces
            [MW0,Mvar0] = calcPowerFromVoltCurrent( ...
                VComplexPhA_busPU,IComplexPhA_busPU,obj.baseMVA);
            assert(all(abs(MW0-obj.genrou_MW0)<obj.absTolCheck))
            assert(all(abs(Mvar0-obj.genrou_Mvar0)<obj.absTolCheck))
        end

        function [obj,vref_out] = excDc1ValInitializations_assumeW0(obj,et_0)
            [x1,x2,x3,x4,x5,vref] = ...
                excDc1ValInitCalculations_assumeW0(obj,et_0,obj.genrou_efd0);
            r = obj.getRowIdx_fromExcDc1Indices(1:obj.numExcDc1)-1;
            obj.x_v(r+1) = x1; obj.x_v(r+2) = x2; obj.x_v(r+3) = x3;
            obj.x_v(r+4) = x4; obj.x_v(r+5) = x5; obj.x_v(r+6) = et_0;
            vref_out = vref;
        end

        function [obj,pref_out] = tgov1dValInitializations_assumeW0(obj)
            % theta = (p/2)*thetaM, so wr = (p/2)*wM. With p = 2 this
            % gives wr = wM, hence Pmech := wM*TorqueM = TorqueM with wr
            % normalised.
            assert(obj.genrouP==2)
            Pmech_TD = obj.genrou_TorqueM0;
            [x1,x2,pref] = tgov1dValInitCalculations_assumeW0(obj,Pmech_TD);
            r = obj.getRowIdx_fromTgov1dIndices(1:obj.numTgov1d)-1;
            obj.x_v(r+1) = x1;
            obj.x_v(r+2) = x2;
            obj.x_v(r+3) = Pmech_TD;
            pref_out = pref;
        end

        function [obj,pset,qset,voDroop_q_set] = ...
                ibrGFMValInitializations_assumeW0(obj,bus_VComplex,bus_Va)
            global idxIbrGFM_IphABC idxIbrGFM_Vg idxIbrGFM_Vo idxIbrGFM_If ...
                   idxIbrGFM_Io idxIbrGFM_Vi idxIbrGFM_Gamma idxIbrGFM_IfRef ...
                   idxIbrGFM_X idxIbrGFM_Wr idxIbrGFM_VoRefDroop_q ...
                   idxIbrGFM_PqAvg idxIbrGFM_Theta idxIbrGFM_Vo_VR_avg ...
                   idxIbrGFM_Vo_VR idxIbrGFM_VoRef

            [vg_dq,vo_dq,vi_dq,if_dq,io_dq,gamma_dq,if_dq_ref,x_dq,wr, ...
                vo_dq_ref,pqavg,iabc_reIm, pset,qset,voDroop_q_set,THETA, ...
                vo_dq_VR_avg,vo_dq_VR,voDroop_q_ref] = ...
                ibrGFMValInitCalculations_assumeW0(obj,bus_VComplex,bus_Va);

            r = obj.getRowIdx_fromIbrGFMIndices(1:obj.numIbrGFM)-1;
            % Interfacing with the network
            for k = 1:6
                obj.x_v(r+idxIbrGFM_IphABC(k)) = iabc_reIm(k,:).';
            end
            % Passive
            obj.x_v(r+idxIbrGFM_Vg(1)) = vg_dq(1,:).';
            obj.x_v(r+idxIbrGFM_Vg(2)) = vg_dq(2,:).';
            obj.x_v(r+idxIbrGFM_Vo(1)) = vo_dq(1,:).';
            obj.x_v(r+idxIbrGFM_Vo(2)) = vo_dq(2,:).';
            obj.x_v(r+idxIbrGFM_If(1)) = if_dq(1,:).';
            obj.x_v(r+idxIbrGFM_If(2)) = if_dq(2,:).';
            obj.x_v(r+idxIbrGFM_Io(1)) = io_dq(1,:).';
            obj.x_v(r+idxIbrGFM_Io(2)) = io_dq(2,:).';
            obj.x_v(r+idxIbrGFM_Vi(1)) = vi_dq(1,:).';
            obj.x_v(r+idxIbrGFM_Vi(2)) = vi_dq(2,:).';
            % Current control
            obj.x_v(r+idxIbrGFM_Gamma(1)) = gamma_dq(1,:).';
            obj.x_v(r+idxIbrGFM_Gamma(2)) = gamma_dq(2,:).';
            % Voltage control
            obj.x_v(r+idxIbrGFM_IfRef(1)) = if_dq_ref(1,:).';
            obj.x_v(r+idxIbrGFM_IfRef(2)) = if_dq_ref(2,:).';
            obj.x_v(r+idxIbrGFM_X(1)) = x_dq(1,:).';
            obj.x_v(r+idxIbrGFM_X(2)) = x_dq(2,:).';
            % Droop
            obj.x_v(r+idxIbrGFM_Wr) = wr;
            obj.x_v(r+idxIbrGFM_VoRefDroop_q) = voDroop_q_ref;
            obj.x_v(r+idxIbrGFM_PqAvg(1)) = pqavg(1,:).';
            obj.x_v(r+idxIbrGFM_PqAvg(2)) = pqavg(2,:).';
            obj.x_v(r+idxIbrGFM_Theta) = THETA;
            % Virtual resistor
            obj.x_v(r+idxIbrGFM_Vo_VR_avg(1)) = vo_dq_VR_avg(1,:);
            obj.x_v(r+idxIbrGFM_Vo_VR_avg(2)) = vo_dq_VR_avg(2,:);
            obj.x_v(r+idxIbrGFM_Vo_VR(1)) = vo_dq_VR(1,:);
            obj.x_v(r+idxIbrGFM_Vo_VR(2)) = vo_dq_VR(2,:);
            obj.x_v(r+idxIbrGFM_VoRef(1)) = vo_dq_ref(1,:);
            obj.x_v(r+idxIbrGFM_VoRef(2)) = vo_dq_ref(2,:);
        end

        function [obj,T0,PQz,PQi,PQp] = ...
                loadCompValInitializations_assumeW0(obj,bus_VComplex)
            global idxLoadComp_Vs idxLoadComp_IphABC idxLoadComp_WSlip ...
                   idxLoadComp_ThetaSys idxLoadComp_Is idxLoadComp_Ir ...
                   idxLoadComp_Wr idxLoadComp_PQ idxLoadComp_Vm idxLoadComp_Izip

            [vdq_s,IComplexPhABC,wslip,thetaSys,wr,idq_s,idq_r, ...
                PQ_Zip,Vm,IComplexPhA_Zip,T0,PQz,PQi,PQp] = ...
                compLoadValInitCalculations_assumeW0(obj,bus_VComplex);

            r = obj.getRowIdx_fromLoadCompIndices(1:obj.numLoadComp)-1;
            obj.x_v(r+idxLoadComp_Vs(1)) = vdq_s(1,:);
            obj.x_v(r+idxLoadComp_Vs(2)) = vdq_s(2,:);
            for k = 1:3
                obj.x_v(r+idxLoadComp_IphABC(2*k-1)) = real(IComplexPhABC(k,:));
                obj.x_v(r+idxLoadComp_IphABC(2*k))   = imag(IComplexPhABC(k,:));
            end
            obj.x_v(r+idxLoadComp_WSlip)    = wslip;
            obj.x_v(r+idxLoadComp_ThetaSys) = thetaSys;
            obj.x_v(r+idxLoadComp_Is(1)) = idq_s(1,:);
            obj.x_v(r+idxLoadComp_Is(2)) = idq_s(2,:);
            obj.x_v(r+idxLoadComp_Ir(1)) = idq_r(1,:);
            obj.x_v(r+idxLoadComp_Ir(2)) = idq_r(2,:);
            obj.x_v(r+idxLoadComp_Wr) = wr;
            obj.x_v(r+idxLoadComp_PQ(1)) = PQ_Zip(1,:);
            obj.x_v(r+idxLoadComp_PQ(2)) = PQ_Zip(2,:);
            obj.x_v(r+idxLoadComp_Vm) = Vm;
            obj.x_v(r+idxLoadComp_Izip(1)) = real(IComplexPhA_Zip);
            obj.x_v(r+idxLoadComp_Izip(2)) = imag(IComplexPhA_Zip);
        end
    end

    methods
        %% ================ Triplet stamp maps (built once) =============
        function obj = buildStampMaps(obj)
            %BUILDSTAMPMAPS  Precompute the row/col index of every triplet.
            %
            % Nothing here depends on x_v, so it runs once. STAMPLINEAR and
            % STAMPNONLINEAR then only refill obj.stampVal, and ASSEMBLE
            % turns the three arrays into the sparse Jacobian.
            global idxGenrou_Idq idxGenrou_Ifd idxGenrou_Edq idxGenrou_Efd ...
                   idxGenrou_Wr idxGenrou_Theta idxGenrou_IphABC idxGenrou_Tm
            global idxExcDc1_x1 idxExcDc1_x2 idxExcDc1_x3 idxExcDc1_x4 ...
                   idxExcDc1_x5 idxExcDc1_Vt
            global idxTgov1d_x1 idxTgov1d_x2 idxTgov1d_Pmech
            global idxIbrGFM_IphABC idxIbrGFM_Vg idxIbrGFM_Vo idxIbrGFM_If ...
                   idxIbrGFM_Io idxIbrGFM_Vi idxIbrGFM_Gamma idxIbrGFM_IfRef ...
                   idxIbrGFM_X idxIbrGFM_Wr idxIbrGFM_VoRefDroop_q ...
                   idxIbrGFM_PqAvg idxIbrGFM_Theta idxIbrGFM_Vo_VR_avg ...
                   idxIbrGFM_Vo_VR idxIbrGFM_VoRef
            global idxLoadComp_Vs idxLoadComp_IphABC idxLoadComp_WSlip ...
                   idxLoadComp_ThetaSys idxLoadComp_Is idxLoadComp_Ir ...
                   idxLoadComp_Wr idxLoadComp_PQ idxLoadComp_Vm idxLoadComp_Izip

            obj = obj.updateTimeSimSamples;

            idxW = obj.getRowIdx_w;
            R = {}; C = {}; slice = struct(); n = 0;
            function addSlice(name,r,c)
                R{end+1} = r(:); C{end+1} = c(:); %#ok<AGROW>
                slice.(name) = n + (1:numel(r)).';
                n = n + numel(r);
            end

            % ---------------- Linear network ---------------------------
            % Transmission lines: 12 phase equations x (12 voltages + w)
            nTx = numel(obj.tx_busFrom);
            [rTx,cTx] = deal(zeros(12*13*nTx,1));
            for k = 1:nTx
                idxFT = obj.branchRowIdxFT(obj.tx_busFrom(k),obj.tx_busTo(k));
                [rr,cc] = saltStampTriplets(idxFT,[idxFT;idxW]);
                rTx((k-1)*156+(1:156)) = rr; cTx((k-1)*156+(1:156)) = cc;
            end
            addSlice('txY',rTx,cTx);

            % Transformers. A lossless transformer carries no admittance,
            % so it gets an explicit branch-current unknown instead (MNA).
            isMna = false(numel(obj.xfmr_busFrom),1);
            isMna(obj.xfmr_mnaIdxMap(:,1)) = true;
            obj.map.xfmr.isMna = isMna;
            obj.map.xfmr.idxAdm = find(~isMna);
            obj.map.xfmr.idxMna = find(isMna);

            rY = []; cY = [];
            for k = obj.map.xfmr.idxAdm.'
                idxFT = obj.branchRowIdxFT(obj.xfmr_busFrom(k),obj.xfmr_busTo(k));
                [rr,cc] = saltStampTriplets(idxFT,[idxFT;idxW]);
                rY = [rY;rr]; cY = [cY;cc]; %#ok<AGROW>
            end
            addSlice('xfmrY',rY,cY);

            rSrc = []; cSrc = []; rZ = []; cZ = [];
            for k = obj.map.xfmr.idxMna.'
                idxFT = obj.branchRowIdxFT(obj.xfmr_busFrom(k),obj.xfmr_busTo(k));
                mna = obj.xfmr_mnaIdxMap(obj.xfmr_mnaIdxMap(:,1)==k,2);
                idxI = obj.getRowIdxAbc_MNA_fromMnaVarIdx(mna);
                % The branch current appears as a source in the two node
                % equations it connects...
                [rr,cc] = saltStampTriplets(idxFT,idxI);
                rSrc = [rSrc;rr]; cSrc = [cSrc;cc]; %#ok<AGROW>
                % ...and gets its own row enforcing V_from - V_to = Z*I
                [rr,cc] = saltStampTriplets(idxI,[idxFT;idxI;idxW]);
                rZ = [rZ;rr]; cZ = [cZ;cc]; %#ok<AGROW>
            end
            addSlice('xfmrIsrc',rSrc,cSrc);
            addSlice('xfmrZ',rZ,cZ);

            % Constant-impedance loads
            rL = []; cL = [];
            for k = 1:numel(obj.load_bus)
                idxF = obj.getRowIdxAbc_fromBusNum(obj.load_bus(k));
                [rr,cc] = saltStampTriplets(idxF,[idxF;idxW]);
                rL = [rL;rr]; cL = [cL;cc]; %#ok<AGROW>
            end
            addSlice('loadY',rL,cL);

            % ---- Current injected by the nonlinear devices -------------
            % Current direction: machines and inverters produce power, so
            % their current flows out of the device into the circuit and is
            % stamped negative. Loads consume power, so their current flows
            % into the device and is stamped positive.
            [rI,cI] = obj.injectionTriplets(obj.genrou_bus, ...
                @(i)obj.getRowIdx_fromGenrouIndices(i),idxGenrou_IphABC);
            addSlice('injGenrou',rI,cI);
            [rI,cI] = obj.injectionTriplets(obj.ibrGFM_bus, ...
                @(i)obj.getRowIdx_fromIbrGFMIndices(i),idxIbrGFM_IphABC);
            addSlice('injIbrGFM',rI,cI);
            [rI,cI] = obj.injectionTriplets(obj.loadComp_bus, ...
                @(i)obj.getRowIdx_fromLoadCompIndices(i),idxLoadComp_IphABC);
            addSlice('injLoadComp',rI,cI);

            % ---------------- GENROU -----------------------------------
            % In-block pattern, in the order the vectorised kernel emits.
            g = idxGenrou_Idq; gf = idxGenrou_Ifd; ge = idxGenrou_Edq;
            gE = idxGenrou_Efd; gw = idxGenrou_Wr; gt = idxGenrou_Theta;
            gi = idxGenrou_IphABC; gT = idxGenrou_Tm;
            P = { 1:3,      [ge gE g gf gw]      % stator
                  4:5,      [gt ge]              % park transform, voltage
                  gw,       gw                   % rotor speed = system speed
                  8,        [gT gw ge g]         % swing
                  9:10,     [gi(1:2) gt g]       % park transform, current
                  11:14,    gi };                % balanced output current
            obj.map.genrou = obj.localPattern(P,obj.numGenrouEqns);
            [r,c] = obj.deviceTriplets(obj.map.genrou, ...
                @(i)obj.getRowIdx_fromGenrouIndices(i),obj.numGenrou);
            % Park-transform voltage couples to the bus, speed to w
            [rp,cp] = obj.busCoupleTriplets(obj.genrou_bus, ...
                @(i)obj.getRowIdx_fromGenrouIndices(i),ge,1:2);
            [rs,cs] = obj.wCoupleTriplets(@(i)obj.getRowIdx_fromGenrouIndices(i), ...
                obj.numGenrou,gw,idxW);
            addSlice('genrou',[r;rp;rs],[c;cp;cs]);

            % ---------------- EXDC1 ------------------------------------
            e1 = idxExcDc1_x1; e2 = idxExcDc1_x2; e3 = idxExcDc1_x3;
            e4 = idxExcDc1_x4; e5 = idxExcDc1_x5; ev = idxExcDc1_Vt;
            P = { e1, [e1 e4]        % -Ke*x1 + x4 = 0
                  e2, [e2 ev]        %  x2 - vt   = 0
                  e3, [e2 e3 e5]     %  x2 + x3 + x5 - vref = 0
                  e4, [e3 e4]        % -Ka*x3 + x4 = 0
                  e5, e5             %  x5 = 0
                  ev, ev };          % -vt + |vabc| = 0
            obj.map.excDc1 = obj.localPattern(P,obj.numExcDc1Eqns);
            [r,c] = obj.deviceTriplets(obj.map.excDc1, ...
                @(i)obj.getRowIdx_fromExcDc1Indices(i),obj.numExcDc1);
            % vt row reads the three real bus-voltage components
            rV = []; cV = []; rF = []; cF = [];
            for k = 1:obj.numExcDc1
                base = obj.getRowIdx_fromExcDc1Indices(k)-1;
                vabc = obj.getRowIdxAbc_fromBusNum(obj.excDc1_bus(k));
                [rr,cc] = saltStampTriplets(base+ev,vabc(1:2:end));
                rV = [rV;rr]; cV = [cV;cc]; %#ok<AGROW>
                % ...and the exciter drives the machine's efd equation
                gbase = obj.getRowIdx_fromGenrouIndices(obj.excDc1_genrouIdx(k))-1;
                rF = [rF; repmat(gbase+gE,3,1)]; %#ok<AGROW>
                cF = [cF; base+e1; gbase+gE; gbase+gw]; %#ok<AGROW>
            end
            addSlice('excDc1',[r;rV;rF],[c;cV;cF]);

            % ---------------- TGOV1D -----------------------------------
            t1 = idxTgov1d_x1; t2 = idxTgov1d_x2; tp = idxTgov1d_Pmech;
            P = { t1, [t1 t2]     % x1 - x2 = 0
                  t2, t2          % pref - (wr/w0-1) - R*x2 = 0
                  tp, [t1 tp] };  % x1 - (wr/w0-1)*Dt - pmech = 0
            obj.map.tgov1d = obj.localPattern(P,obj.numTgov1dEqns);
            [r,c] = obj.deviceTriplets(obj.map.tgov1d, ...
                @(i)obj.getRowIdx_fromTgov1dIndices(i),obj.numTgov1d);
            rW = []; cW = []; rT = []; cT = [];
            for k = 1:obj.numTgov1d
                base  = obj.getRowIdx_fromTgov1dIndices(k)-1;
                gbase = obj.getRowIdx_fromGenrouIndices(obj.tgov1d_genrouIdx(k))-1;
                % All three governor rows see the machine speed
                rW = [rW; base+(1:3).']; cW = [cW; repmat(gbase+gw,3,1)]; %#ok<AGROW>
                % ...and the governor drives the machine's torque equation
                rT = [rT; repmat(gbase+gT,3,1)]; %#ok<AGROW>
                cT = [cT; base+tp; gbase+gw; gbase+gT]; %#ok<AGROW>
            end
            addSlice('tgov1d',[r;rW;rT],[c;cW;cT]);

            % ---------------- GFM IBR ----------------------------------
            bI = idxIbrGFM_IphABC; bVg = idxIbrGFM_Vg; bVo = idxIbrGFM_Vo;
            bIf = idxIbrGFM_If; bIo = idxIbrGFM_Io; bVi = idxIbrGFM_Vi;
            bG = idxIbrGFM_Gamma; bIr = idxIbrGFM_IfRef; bX = idxIbrGFM_X;
            bw = idxIbrGFM_Wr; bq = idxIbrGFM_VoRefDroop_q;
            bpq = idxIbrGFM_PqAvg; bth = idxIbrGFM_Theta;
            bva = idxIbrGFM_Vo_VR_avg; bvr = idxIbrGFM_Vo_VR;
            bvR = idxIbrGFM_VoRef;
            P = { 1:2,   [bI(1:2) bth bIo]              % park transform, current
                  3:6,   bI                             % balanced output current
                  7:8,   [bth bVg]                      % park transform, voltage
                  9:10,  [bIo bw bVo bVg]               % passive coupling
                  11:12, [bIf bIo bVo bw]               % passive capacitor
                  13:14, [bIf bw bVi bVo]               % passive filter
                  15:16, [bVi bG bIf bw bVo]            % current controller
                  17:18, [bIr bIf]                      % current integral error
                  19:20, [bIr bX bVo bw bIo]            % voltage controller
                  21:22, [bvR bVo]                      % voltage integral error
                  bw,    [bw bpq(1)]                    % P droop
                  bq,    [bq bpq(2)]                    % Q droop
                  bpq(1),[bpq(1) bVg bIo]               % measured P
                  bpq(2),[bpq(2) bVg bIo]               % measured Q
                  bth,   bw                             % w_sys = wr
                  bva,   [bva bIo]                      % virtual resistor, filtered
                  bvr,   [bvr bva bIo]                  % virtual resistor drop
                  bvR,   [bvR bvr]                      % output voltage reference
                  bvR(2),bq };                          % ...plus the Q droop term
            obj.map.ibrGFM = obj.localPattern(P,obj.numIbrGFMEqns);
            [r,c] = obj.deviceTriplets(obj.map.ibrGFM, ...
                @(i)obj.getRowIdx_fromIbrGFMIndices(i),obj.numIbrGFM);
            [rp,cp] = obj.busCoupleTriplets(obj.ibrGFM_bus, ...
                @(i)obj.getRowIdx_fromIbrGFMIndices(i),bVg,1:2);
            [rs,cs] = obj.wCoupleTriplets(@(i)obj.getRowIdx_fromIbrGFMIndices(i), ...
                obj.numIbrGFM,bth,idxW);
            addSlice('ibrGFM',[r;rp;rs],[c;cp;cs]);

            % ---------------- Composite load ---------------------------
            lv = idxLoadComp_Vs; li = idxLoadComp_IphABC;
            ls = idxLoadComp_WSlip; lt = idxLoadComp_ThetaSys;
            lIs = idxLoadComp_Is; lIr = idxLoadComp_Ir; lw = idxLoadComp_Wr;
            lpq = idxLoadComp_PQ; lvm = idxLoadComp_Vm; lz = idxLoadComp_Izip;
            P = { 1:2,   lv                      % park transform, voltage
                  3:4,   [lIs lz li(1:2)]        % park transform, current
                  5:8,   li                      % balanced output current
                  ls,    [ls lw]                 % slip frequency
                  lt,    lt                      % system angle
                  11:14, [lv lIs lIr ls]         % stator / rotor
                  lw,    [lw lIs lIr]            % swing
                  16:17, [lpq lvm]               % ZIP power
                  lvm,   [lvm lv]                % ZIP voltage magnitude
                  19:20, [lpq lv lvm lz] };      % ZIP current
            obj.map.loadComp = obj.localPattern(P,obj.numLoadCompEqns);
            [r,c] = obj.deviceTriplets(obj.map.loadComp, ...
                @(i)obj.getRowIdx_fromLoadCompIndices(i),obj.numLoadComp);
            [rp,cp] = obj.busCoupleTriplets(obj.loadComp_bus, ...
                @(i)obj.getRowIdx_fromLoadCompIndices(i),lv,1:2);
            % The slip, stator and swing equations all see the system speed
            [r1,c1] = obj.wCoupleTriplets(@(i)obj.getRowIdx_fromLoadCompIndices(i), ...
                obj.numLoadComp,ls,idxW);
            [r2,c2] = obj.wCoupleTriplets(@(i)obj.getRowIdx_fromLoadCompIndices(i), ...
                obj.numLoadComp,lIs,idxW);
            [r3,c3] = obj.wCoupleTriplets(@(i)obj.getRowIdx_fromLoadCompIndices(i), ...
                obj.numLoadComp,lw,idxW);
            addSlice('loadComp',[r;rp;r1;r2;r3],[c;cp;c1;c2;c3]);

            % ---------------- Slack-angle (frequency) equation ----------
            vSlack = obj.getRowIdxAbc_fromBusNum(obj.slack_busNum);
            addSlice('phase',[idxW;idxW],vSlack(1:2));

            % ---- Row maps for the right-hand side ----------------------
            % The RHS is accumulated, so each branch just needs the list of
            % equations its 12 (or 6) real entries land on.
            rowsTx = zeros(12*nTx,1);
            for k = 1:nTx
                rowsTx((k-1)*12+(1:12)) = ...
                    obj.branchRowIdxFT(obj.tx_busFrom(k),obj.tx_busTo(k));
            end
            obj.map.txIrows = rowsTx;
            obj.map.xfmrYIrows = cell2mat(arrayfun(@(k) ...
                obj.branchRowIdxFT(obj.xfmr_busFrom(k),obj.xfmr_busTo(k)), ...
                obj.map.xfmr.idxAdm,'UniformOutput',false));
            obj.map.xfmrZIrows = cell2mat(arrayfun(@(k) ...
                obj.getRowIdxAbc_MNA_fromMnaVarIdx( ...
                obj.xfmr_mnaIdxMap(obj.xfmr_mnaIdxMap(:,1)==k,2)), ...
                obj.map.xfmr.idxMna,'UniformOutput',false));
            obj.map.loadIrows = cell2mat(arrayfun(@(b) ...
                obj.getRowIdxAbc_fromBusNum(b), ...
                obj.load_bus(:),'UniformOutput',false));

            % ---------------- Finalise ---------------------------------
            obj.stampRow = vertcat(R{:});
            obj.stampCol = vertcat(C{:});
            obj.stampVal = zeros(n,1);
            obj.stampSlice = slice;

            % Entries that never change are filled once here
            obj.stampVal(slice.xfmrIsrc) = ...
                repmat(reshape(blkdiag([eye(2);-eye(2)],[eye(2);-eye(2)], ...
                [eye(2);-eye(2)]),[],1),numel(obj.map.xfmr.idxMna),1);
            obj.stampVal(slice.injGenrou)   = repmat(reshape(-eye(6),[],1),obj.numGenrou,1);
            obj.stampVal(slice.injIbrGFM)   = repmat(reshape(-eye(6),[],1),obj.numIbrGFM,1);
            obj.stampVal(slice.injLoadComp) = repmat(reshape( eye(6),[],1),obj.numLoadComp,1);
            obj.stampVal(slice.phase)       = [obj.slack_imScalConst; -1];
        end
    end

    methods (Hidden)
        %% ------------- Small helpers used by BUILDSTAMPMAPS -----------
        function idxFT = branchRowIdxFT(obj,busFrom,busTo)
            % Interleave the from/to node rows phase by phase, so a branch
            % block reads [Fa_re Fa_im Ta_re Ta_im Fb_re ... ].
            idxF = obj.getRowIdxAbc_fromBusNum(busFrom);
            idxT = obj.getRowIdxAbc_fromBusNum(busTo);
            idxFT = reshape([idxF(1:2:end),idxF(2:2:end), ...
                             idxT(1:2:end),idxT(2:2:end)].',[],1);
        end

        function m = localPattern(~,P,nEqns)
            % Turn a {rows,cols;...} block list into flat local row/col
            % index vectors plus the selector that sums value contributions
            % back onto equations.
            row = []; col = [];
            for k = 1:size(P,1)
                [rr,cc] = saltStampTriplets(P{k,1},P{k,2});
                row = [row;rr]; col = [col;cc]; %#ok<AGROW>
            end
            m.row = row; m.col = col; m.nVals = numel(row);
            m.S = sparse(row,1:numel(row),1,nEqns,numel(row));
        end

        function [r,c] = deviceTriplets(~,m,rowBaseFcn,nDev)
            % Replicate one device's local pattern across every device.
            if nDev==0, r = []; c = []; return, end
            base = rowBaseFcn((1:nDev).') - 1;          % [nDev x 1]
            r = reshape(m.row + base.',[],1);
            c = reshape(m.col + base.',[],1);
        end

        function [r,c] = busCoupleTriplets(obj,busList,rowBaseFcn,locRows,phaseCols)
            % Device rows that read the terminal bus voltage directly.
            r = []; c = [];
            for k = 1:numel(busList)
                base = rowBaseFcn(k)-1;
                vabc = obj.getRowIdxAbc_fromBusNum(busList(k));
                [rr,cc] = saltStampTriplets(base+locRows,vabc(phaseCols));
                r = [r;rr]; c = [c;cc]; %#ok<AGROW>
            end
        end

        function [r,c] = wCoupleTriplets(~,rowBaseFcn,nDev,locRows,idxW)
            % Device rows that read the system frequency.
            if nDev==0, r = []; c = []; return, end
            base = rowBaseFcn((1:nDev).')-1;
            r = reshape(base.' + locRows(:),[],1);
            c = repmat(idxW,numel(r),1);
        end

        function [r,c] = injectionTriplets(obj,busList,rowBaseFcn,locI)
            % Bus equations that receive a device's terminal current.
            r = []; c = [];
            for k = 1:numel(busList)
                vabc = obj.getRowIdxAbc_fromBusNum(busList(k));
                iabc = rowBaseFcn(k)-1+locI;
                [rr,cc] = saltStampTriplets(vabc,iabc);
                r = [r;rr]; c = [c;cc]; %#ok<AGROW>
            end
        end
    end

    methods
        %% ====================== Stamping ==============================
        function obj = stampLinear_Y_I(obj)
            %STAMPLINEAR_Y_I  Values for the passive network, all at once.
            sysSize = obj.getSysSize;
            idxW = obj.getRowIdx_w;
            obj.I = zeros(sysSize,1);

            w = obj.x_v(idxW);

            % ---- Transmission lines (3ph PI, no mutual coupling) -------
            % Eqn at the "from" node:  Ifout + (Vf-Vt)*YeqBr + Vf*YeqSh = 0
            % Eqn at the "to"   node:  Itout + (Vt-Vf)*YeqBr + Vt*YeqSh = 0
            nTx = numel(obj.tx_busFrom);
            if nTx > 0
                [VF,VT] = obj.branchVoltages(obj.tx_busFrom,obj.tx_busTo);
                Yb  = (1./(obj.tx_R + 1j*w*obj.tx_L)).';       % 1 x nTx
                Ysh = (1j*w*obj.tx_C).';
                Ibr = Yb.*(VF-VT);
                Iabc = [Ibr + Ysh.*VF; -Ibr + Ysh.*VT];
                % dI/dV: series branch plus the shunt on each end
                dIdV = SALT.pairBlock2x2(Yb+Ysh,-Yb,-Yb,Yb+Ysh);
                % dI/dw: series branch plus both shunts
                dBr = SALT.calc_dIdw_RLSeries_vec(obj.tx_R.',obj.tx_L.',VF-VT,w);
                dIdw = [dBr + 1j*obj.tx_C.'.*VF; -dBr + 1j*obj.tx_C.'.*VT];
                [obj.stampVal(obj.stampSlice.txY),Ivals] = ...
                    SALT.networkBlock(dIdV,SALT.interleaveFT(dIdw), ...
                    SALT.interleaveFT([VF;VT]),SALT.interleaveFT(Iabc),w);
                obj.I = obj.accumI(obj.I,obj.map.txIrows,Ivals);
            end

            % ---- Transformers with a nonzero resistance (admittance) ---
            iA = obj.map.xfmr.idxAdm;
            if ~isempty(iA)
                [VF,VT] = obj.branchVoltages(obj.xfmr_busFrom(iA),obj.xfmr_busTo(iA));
                Ym = (1./(obj.xfmr_R(iA) + 1j*w*obj.xfmr_L(iA))).';
                Iabc3 = Ym.*(VF-VT);
                dIdV = SALT.pairBlock2x2(Ym,-Ym,-Ym,Ym);
                dIdw = SALT.calc_dIdw_RLSeries_vec(obj.xfmr_R(iA).', ...
                    obj.xfmr_L(iA).',VF-VT,w);
                [obj.stampVal(obj.stampSlice.xfmrY),Ivals] = ...
                    SALT.networkBlock(dIdV, ...
                    SALT.interleaveFT([dIdw;-dIdw]), ...
                    SALT.interleaveFT([VF;VT]), ...
                    SALT.interleaveFT([Iabc3;-Iabc3]),w);
                obj.I = obj.accumI(obj.I,obj.map.xfmrYIrows,Ivals);
            end

            % ---- Lossless transformers (branch current as an unknown) --
            iM = obj.map.xfmr.idxMna;
            if ~isempty(iM)
                [VF,VT] = obj.branchVoltages(obj.xfmr_busFrom(iM),obj.xfmr_busTo(iM));
                Iabc = obj.mnaCurrents(iM);
                nM = numel(iM);
                L = obj.xfmr_L(iM).';
                Z = 1j*w*L;                                    % 1 x nM
                % F = -(Vf - Vt) + Z*I
                F = -(VF-VT) + Z.*Iabc;
                dFdV = repmat(-blkdiag([1 -1],[1 -1],[1 -1]),1,1,nM);
                dFdI = reshape(eye(3),3,3,1).*reshape(Z,1,1,nM);
                % dZ/dw = jL, so dF/dw = jL*I
                dFdw = reshape(1j*L.*Iabc,3,1,nM);
                M = cat(2,dFdV,dFdI,dFdw);
                E = SALT.reImDropLast(M);
                obj.stampVal(obj.stampSlice.xfmrZ) = E(:);
                Ic = pagemtimes(dFdV,reshape(SALT.interleaveFT([VF;VT]),6,1,nM)) ...
                   + pagemtimes(dFdI,reshape(Iabc,3,1,nM)) ...
                   + dFdw*w - reshape(F,3,1,nM);
                obj.I = obj.accumI(obj.I,obj.map.xfmrZIrows, ...
                    SALT.interleaveReIm(reshape(Ic,3,nM)));
            end

            % ---- Constant-impedance loads ------------------------------
            nLd = numel(obj.load_bus);
            if nLd > 0
                idx = obj.getRowIdxAbc_fromBusNums(obj.load_bus);
                % getRowIdxAbc_fromBusNums returns phase-major
                % ([phA(all buses); phB(all); phC(all)]), so the gathered
                % values run bus-fastest within each phase: reshape to
                % (nLd x 3) and transpose, NOT straight to (3 x nLd).
                Vabc = reshape(complex(obj.x_v(idx(1:2:end)), ...
                    obj.x_v(idx(2:2:end))),nLd,3).';
                [Ym,dYdw] = obj.constLoadAdmittance(w);
                Iabc = Ym.*Vabc;
                dIdV = reshape(eye(3),3,3,1).*reshape(Ym,1,1,nLd);
                dIdw = reshape(dYdw.*Vabc,3,1,nLd);
                M = cat(2,dIdV,dIdw);
                E = SALT.reImDropLast(M);
                obj.stampVal(obj.stampSlice.loadY) = E(:);
                Ic = pagemtimes(dIdV,reshape(Vabc,3,1,nLd)) + dIdw*w ...
                   - reshape(Iabc,3,1,nLd);
                obj.I = obj.accumI(obj.I,obj.map.loadIrows, ...
                    SALT.interleaveReIm(reshape(Ic,3,nLd)));
            end

        end

        function stampPhase_Y_I(obj)
            % Extra equation that pins the slack angle, needed because the
            % system frequency is itself an unknown:
            %   0 = (CONST)*V_real + (-1)*V_imag
            obj.I(obj.getRowIdx_w) = 0;
        end

        function obj = assemble_Y_I(obj)
            obj.Y = saltAssembleY(obj.stampRow,obj.stampCol,obj.stampVal, ...
                obj.getSysSize);
        end
    end

    methods (Hidden)
        function [VF,VT] = branchVoltages(obj,busFrom,busTo)
            n = numel(busFrom);
            iF = obj.getRowIdxAbc_fromBusNums(busFrom);
            iT = obj.getRowIdxAbc_fromBusNums(busTo);
            % getRowIdxAbc_fromBusNums returns [phA(all); phB(all); phC(all)]
            VF = reshape(complex(obj.x_v(iF(1:2:end)),obj.x_v(iF(2:2:end))),n,3).';
            VT = reshape(complex(obj.x_v(iT(1:2:end)),obj.x_v(iT(2:2:end))),n,3).';
        end

        function Iabc = mnaCurrents(obj,deviceIdx)
            n = numel(deviceIdx);
            Iabc = zeros(3,n);
            for k = 1:n
                mna = obj.xfmr_mnaIdxMap(obj.xfmr_mnaIdxMap(:,1)==deviceIdx(k),2);
                idxI = obj.getRowIdxAbc_MNA_fromMnaVarIdx(mna);
                Iabc(:,k) = complex(obj.x_v(idxI(1:2:end)),obj.x_v(idxI(2:2:end)));
            end
        end

        function [Ym,dYdw] = constLoadAdmittance(obj,w)
            n = numel(obj.load_bus);
            Ym = zeros(1,n); dYdw = zeros(1,n);
            for k = 1:n
                X = obj.load_X(k); Rr = obj.load_R(k);
                if X > 0
                    L = X/obj.w0;
                    Ym(k) = 1/(Rr + 1j*w*L);
                    dYdw(k) = -1j*L/(Rr + 1j*w*L)^2;
                elseif X < 0
                    Rp = 1/real(1/obj.load_Z(k));
                    Cc = imag(1/obj.load_Z(k))/obj.w0;
                    Ym(k) = 1/Rp + 1j*w*Cc;
                    dYdw(k) = 1j*Cc;
                else
                    Ym(k) = 1/Rr;
                    dYdw(k) = 0;
                end
            end
        end

        function Iout = accumI(~,Iout,rows,vals)
            Iout = Iout + accumarray(rows(:),vals(:),size(Iout));
        end
    end

    methods (Static, Hidden)
        function B = pairBlock2x2(a11,a12,a21,a22)
            % Build the 6x6xN block-diagonal [a11 a12; a21 a22] repeated
            % once per phase, from N-length rows of scalars.
            n = numel(a11);
            B = zeros(6,6,n);
            for p = 0:2
                B(2*p+1,2*p+1,:) = a11;  B(2*p+1,2*p+2,:) = a12;
                B(2*p+2,2*p+1,:) = a21;  B(2*p+2,2*p+2,:) = a22;
            end
        end

        function X = interleaveFT(X6)
            % [F(3); T(3)] per column -> [Fa;Ta;Fb;Tb;Fc;Tc]
            X = X6([1 4 2 5 3 6],:);
        end

        function v = interleaveReIm(Xc)
            % p x n complex -> 2p x n real with (re,im) interleaved, then
            % flattened column-major.
            [p,n] = size(Xc);
            v = zeros(2*p,n);
            v(1:2:end,:) = real(Xc);
            v(2:2:end,:) = imag(Xc);
            v = v(:);
        end

        function E = reImDropLast(M)
            % p x q x n complex -> 2p x (2q-1) x n real. The final column
            % is dropped because the extra variable (system frequency) is
            % real, so only its real-part column is meaningful.
            [p,q,n] = size(M);
            E = zeros(2*p,2*q,n);
            E(1:2:end,1:2:end,:) =  real(M);
            E(1:2:end,2:2:end,:) = -imag(M);
            E(2:2:end,1:2:end,:) =  imag(M);
            E(2:2:end,2:2:end,:) =  real(M);
            E(:,end,:) = [];
        end

        function [vals,Ivals] = networkBlock(dIdV,dIdw,Vstamp,Istamp,w)
            % Shared assembly for every passive branch: expand
            % [dI/dV, dI/dw] into the real-imaginary stamp and form the
            % Newton right-hand side J*x - F.
            n = size(dIdV,3);
            p = size(dIdV,1);
            M = cat(2,dIdV,reshape(dIdw,p,1,n));
            E = SALT.reImDropLast(M);
            vals = E(:);
            Ic = pagemtimes(dIdV,reshape(Vstamp,p,1,n)) ...
               + reshape(dIdw,p,1,n)*w - reshape(Istamp,p,1,n);
            Ivals = SALT.interleaveReIm(reshape(Ic,p,n));
        end

        function dI_dw = calc_dIdw_RLSeries_vec(R,L,V,w)
            % d/dw of I = V/(R + jwL), evaluated per column of V.
            absZsq = (R.^2 + (w*L).^2);
            dI_dw = -( (1j*L)./absZsq + (2*w*L.^2).*(R-1j*w*L)./absZsq.^2 ).*V;
        end

        function Yeq = calcYeq_RLSeries(R,L,w),  Yeq = 1/(R + 1j*w*L); end
        function Yeq = calcYeq_C(C,w),           Yeq = 1j*w*C;         end
        function Zeq = calcZeq_L(L,w),           Zeq = 1j*w*L;         end
    end

    methods
        %% ============== Nonlinear device stamping =====================
        function obj = stampNonlinear_Y_I(obj)
            %STAMPNONLINEAR_Y_I  One vectorised pass per device model.
            obj = obj.stampGenrou;
            obj = obj.stampExcDc1;
            obj = obj.stampTgov1d;
            obj = obj.stampIbrGFM;
            obj = obj.stampLoadComp;
        end
    end

    methods (Hidden)
        function P = gatherPages(obj,rowBase,loc)
            % x_v entries for one local index set, as [nLoc x 1 x nDev]
            X = obj.x_v(rowBase(:).' + loc(:));       % nLoc x nDev
            P = reshape(X,numel(loc),1,numel(rowBase));
        end

        function V = terminalVphA(obj,busList)
            % Phase-A (re,im) bus voltage of each device, as [2 x 1 x nDev]
            n = numel(busList);
            idx = obj.getRowIdxAbc_fromBusNums(busList);
            V = reshape(obj.x_v(idx(1:2*n)),2,1,n);
        end

        function obj = writeBlockRhs(obj,rowBase,nEqns,vals)
            rows = reshape(rowBase(:).' + (1:nEqns).',[],1);
            obj.I(rows) = vals(:);
        end
        function obj = addBlockRhs(obj,rowBase,locRows,vals)
            rows = reshape(rowBase(:).' + locRows(:),[],1);
            obj.I(rows) = obj.I(rows) + vals(:);
        end

        function Jx = blockJx(obj,m,rowBase,Jvals)
            % J(x)*x for every device, straight from the triplets.
            xLoc = obj.x_v(m.col + rowBase(:).');     % nVals x nDev
            Jx = m.S * (Jvals .* xLoc);               % nEqns x nDev
        end

        %% ---------------------------- GENROU -------------------------
        function obj = stampGenrou(obj)
            n = obj.numGenrou;
            if n==0, return, end
            global idxGenrou_Idq idxGenrou_Ifd idxGenrou_Edq idxGenrou_Efd ...
                   idxGenrou_Wr idxGenrou_Theta idxGenrou_IphABC idxGenrou_Tm

            r  = obj.getRowIdx_fromGenrouIndices((1:n).')-1;
            ws = obj.x_v(obj.getRowIdx_w);
            [~,busIdx] = ismember(obj.genrou_bus,obj.bus);
            toBusI = reshape(obj.genrouToBus_I(busIdx,(1:n).'),1,1,n);
            toBusV = reshape(obj.genrouToBus_V(busIdx,(1:n).')*ones(n,1),1,1,n);

            [Jvals,F,dParkV,dFspeed_dWs] = stampOfGenrouMatrices_SS( ...
                obj.gatherPages(r,[idxGenrou_Idq idxGenrou_Ifd]), ...
                obj.gatherPages(r,[idxGenrou_Edq idxGenrou_Efd]), ...
                obj.gatherPages(r,idxGenrou_IphABC), ...
                obj.terminalVphA(obj.genrou_bus), ...
                obj.gatherPages(r,idxGenrou_Theta), ...
                obj.gatherPages(r,idxGenrou_Wr), ...
                obj.gatherPages(r,idxGenrou_Tm), ws, ...
                obj.genrou_RMatrix,obj.genrou_LspMatrix, ...
                obj.w0, toBusI,toBusV);

            obj.stampVal(obj.stampSlice.genrou) = [Jvals(:);dParkV(:);dFspeed_dWs(:)];

            % J(x[k])*x[k+1] = J(x[k])*x[k] - F(x[k])
            obj = obj.writeBlockRhs(r,obj.numGenrouEqns, ...
                obj.blockJx(obj.map.genrou,r,Jvals) - F);
            % Park transform reads the terminal voltage
            VphA = obj.terminalVphA(obj.genrou_bus);
            obj = obj.addBlockRhs(r,idxGenrou_Edq, ...
                reshape(pagemtimes(reshape(dParkV,2,2,n),VphA),2,n));
            % Speed equation reads the system frequency
            obj = obj.addBlockRhs(r,idxGenrou_Wr,dFspeed_dWs*ws);
        end

        %% ---------------------------- EXDC1 --------------------------
        function obj = stampExcDc1(obj)
            n = obj.numExcDc1;
            if n==0, return, end
            global idxGenrou_Efd idxGenrou_Wr idxExcDc1_x1 idxExcDc1_Vt

            r  = obj.getRowIdx_fromExcDc1Indices((1:n).')-1;
            rg = obj.getRowIdx_fromGenrouIndices(obj.excDc1_genrouIdx(:))-1;
            [~,busIdx] = ismember(obj.excDc1_bus,obj.bus);
            toGenV = obj.busToGenrou_V(busIdx,(1:n).')*ones(n,1);

            % Real parts of the three terminal-voltage phases
            idxVabc = obj.getRowIdxAbc_fromBusNums(obj.excDc1_bus);
            vabcRe = reshape(obj.x_v(idxVabc(1:2:end)),n,3).';   % 3 x n

            [Jvals,Ilin, Fvt,Jvt_dvt,Jvt_dvabc, ...
                Fefd,Jefd_x1,Jefd_efd,Jefd_wr] = stampOfExcDc1Matrices_SS( ...
                obj.x_v(r+idxExcDc1_x1), obj.x_v(rg+idxGenrou_Efd), ...
                obj.x_v(r+idxExcDc1_Vt), vabcRe, obj.x_v(rg+idxGenrou_Wr), ...
                obj.excDc1_Ka,obj.excDc1_Ke,obj.excDc1_vref, ...
                obj.w0, toGenV,obj.excDc1_genrouToExcDc1);

            obj.stampVal(obj.stampSlice.excDc1) = ...
                [reshape([Jvals;Jvt_dvt],[],1); Jvt_dvabc(:); ...
                 reshape([Jefd_x1;Jefd_efd;Jefd_wr],[],1)];

            % Linear state rows carry only the constant right-hand side
            rowsLin = reshape(r(:).' + (1:5).',[],1);
            obj.I(rowsLin) = Ilin(:);
            % Terminal-voltage row
            % J(x[k])*x[k+1] = J(x[k])*x[k] - F(x[k])
            obj.I(r+idxExcDc1_Vt) = (Jvt_dvt.*obj.x_v(r+idxExcDc1_Vt).' ...
                + sum(Jvt_dvabc.*vabcRe,1) - Fvt).';
            % Field-voltage row of the machine this exciter drives
            obj.I(rg+idxGenrou_Efd) = ( ...
                  Jefd_x1 .*obj.x_v(r+idxExcDc1_x1).' ...
                + Jefd_efd.*obj.x_v(rg+idxGenrou_Efd).' ...
                + Jefd_wr .*obj.x_v(rg+idxGenrou_Wr).' - Fefd).';
        end

        %% ---------------------------- TGOV1D -------------------------
        function obj = stampTgov1d(obj)
            n = obj.numTgov1d;
            if n==0, return, end
            global idxGenrou_Wr idxGenrou_Tm idxTgov1d_Pmech

            r  = obj.getRowIdx_fromTgov1dIndices((1:n).')-1;
            rg = obj.getRowIdx_fromGenrouIndices(obj.tgov1d_genrouIdx(:))-1;
            pmech = obj.x_v(r+idxTgov1d_Pmech);
            Tm    = obj.x_v(rg+idxGenrou_Tm);
            wr    = obj.x_v(rg+idxGenrou_Wr);

            [Jvals,Jwr,Ilin,FTm,JTm] = stampOfTgov1dMatrices_SS( ...
                pmech,Tm,wr, obj.tgov1d_R,obj.tgov1d_Dt,obj.tgov1d_Pref, ...
                obj.w0, obj.genrouP);

            obj.stampVal(obj.stampSlice.tgov1d) = [Jvals(:);Jwr(:);JTm(:)];

            obj = obj.writeBlockRhs(r,obj.numTgov1dEqns,Ilin);
            % Mechanical-torque row of the machine this governor drives
            % J(x[k])*x[k+1] = J(x[k])*x[k] - F(x[k])
            obj.I(rg+idxGenrou_Tm) = ...
                (sum(JTm.*[pmech.';wr.';Tm.'],1) - FTm).';
        end

        %% --------------------------- GFM IBR --------------------------
        function obj = stampIbrGFM(obj)
            n = obj.numIbrGFM;
            if n==0, return, end
            global idxIbrGFM_IphABC idxIbrGFM_Vg idxIbrGFM_Vo idxIbrGFM_If ...
                   idxIbrGFM_Io idxIbrGFM_Vi idxIbrGFM_Gamma idxIbrGFM_IfRef ...
                   idxIbrGFM_X idxIbrGFM_Wr idxIbrGFM_VoRefDroop_q ...
                   idxIbrGFM_PqAvg idxIbrGFM_Theta idxIbrGFM_Vo_VR_avg ...
                   idxIbrGFM_Vo_VR idxIbrGFM_VoRef

            r  = obj.getRowIdx_fromIbrGFMIndices((1:n).')-1;
            ws = obj.x_v(obj.getRowIdx_w);
            [~,busIdx] = ismember(obj.ibrGFM_bus,obj.bus);
            toBusI = reshape(obj.ibrGFMToBus_I(busIdx,(1:n).'),1,1,n);
            toBusV = reshape(obj.ibrGFMToBus_V(busIdx,(1:n).'),1,1,n);
            p3 = @(v) reshape(v,1,1,n);

            [Jvals,F,dParkV,dFspeed_dWs] = stampOfIbrGFMMatrices_SS( ...
                obj.gatherPages(r,idxIbrGFM_IphABC), ...
                obj.gatherPages(r,idxIbrGFM_Vg), obj.gatherPages(r,idxIbrGFM_Vo), ...
                obj.gatherPages(r,idxIbrGFM_If), obj.gatherPages(r,idxIbrGFM_Io), ...
                obj.gatherPages(r,idxIbrGFM_Vi), obj.gatherPages(r,idxIbrGFM_Gamma), ...
                obj.gatherPages(r,idxIbrGFM_IfRef), obj.gatherPages(r,idxIbrGFM_X), ...
                obj.gatherPages(r,idxIbrGFM_Wr), obj.gatherPages(r,idxIbrGFM_VoRef), ...
                obj.gatherPages(r,idxIbrGFM_PqAvg(1)), ...
                obj.gatherPages(r,idxIbrGFM_PqAvg(2)), ...
                obj.gatherPages(r,idxIbrGFM_Theta), ...
                obj.terminalVphA(obj.ibrGFM_bus), ...
                obj.gatherPages(r,idxIbrGFM_Vo_VR_avg), ...
                obj.gatherPages(r,idxIbrGFM_Vo_VR), ...
                obj.gatherPages(r,idxIbrGFM_VoRefDroop_q), ws, ...
                p3(obj.ibrGFM_Lf),p3(obj.ibrGFM_Rf),p3(obj.ibrGFM_Cf), ...
                p3(obj.ibrGFM_Lc),p3(obj.ibrGFM_Rc),p3(obj.ibrGFM_kC_i), ...
                p3(obj.ibrGFM_GC),p3(obj.ibrGFM_kV_i),p3(obj.ibrGFM_GV), ...
                p3(obj.ibrGFM_Mp),p3(obj.ibrGFM_Mq), ...
                p3(obj.ibrGFM_pset),p3(obj.ibrGFM_qset), ...
                p3(obj.ibrGFM_voDroop_q_set), ...
                p3(obj.ibrGFM_w_VR),p3(obj.ibrGFM_Rv), ...
                obj.w0, toBusI,toBusV);

            obj.stampVal(obj.stampSlice.ibrGFM) = ...
                [Jvals(:);dParkV(:);dFspeed_dWs(:)];

            obj = obj.writeBlockRhs(r,obj.numIbrGFMEqns, ...
                obj.blockJx(obj.map.ibrGFM,r,Jvals) - F);
            VphA = obj.terminalVphA(obj.ibrGFM_bus);
            obj = obj.addBlockRhs(r,idxIbrGFM_Vg, ...
                reshape(pagemtimes(reshape(dParkV,2,2,n),VphA),2,n));
            obj = obj.addBlockRhs(r,idxIbrGFM_Theta,dFspeed_dWs*ws);
        end

        %% ------------------------ Composite load ----------------------
        function obj = stampLoadComp(obj)
            n = obj.numLoadComp;
            if n==0, return, end
            global idxLoadComp_Vs idxLoadComp_IphABC idxLoadComp_WSlip ...
                   idxLoadComp_ThetaSys idxLoadComp_Is idxLoadComp_Ir ...
                   idxLoadComp_Wr idxLoadComp_PQ idxLoadComp_Vm idxLoadComp_Izip

            r  = obj.getRowIdx_fromLoadCompIndices((1:n).')-1;
            ws = obj.x_v(obj.getRowIdx_w);
            [~,busIdx] = ismember(obj.loadComp_bus,obj.bus);
            toI = reshape(obj.busToLoadComp_I(busIdx,(1:n).'),1,1,n);
            toV = reshape(obj.busToLoadComp_V(busIdx,(1:n).')*ones(n,1),1,1,n);
            p3 = @(v) reshape(v,1,1,n);

            [Jvals,F,dParkV,dFslip_dwSys,dFs_dwSys,dFswing_dwSys] = ...
                stampOfCompositeLoadMatrices_SS( ...
                obj.gatherPages(r,idxLoadComp_Vs), ...
                obj.gatherPages(r,idxLoadComp_IphABC), ...
                obj.gatherPages(r,idxLoadComp_WSlip), ...
                obj.gatherPages(r,idxLoadComp_ThetaSys), ...
                obj.gatherPages(r,idxLoadComp_Is), ...
                obj.gatherPages(r,idxLoadComp_Ir), ...
                obj.gatherPages(r,idxLoadComp_Wr), ...
                obj.terminalVphA(obj.loadComp_bus), ws, ...
                obj.gatherPages(r,idxLoadComp_PQ), ...
                obj.gatherPages(r,idxLoadComp_Vm), ...
                obj.gatherPages(r,idxLoadComp_Izip), ...
                obj.loadCompIM_R,obj.loadCompIM_L,obj.loadCompIM_LTe, ...
                p3(obj.loadCompIM_Tm0),p3(obj.loadCompIM_m), ...
                reshape(obj.loadCompZip_PQz,2,1,n), ...
                reshape(obj.loadCompZip_PQi,2,1,n), ...
                reshape(obj.loadCompZip_PQp,2,1,n), ...
                obj.w0, toI,toV);

            obj.stampVal(obj.stampSlice.loadComp) = [Jvals(:);dParkV(:); ...
                dFslip_dwSys(:);dFs_dwSys(:);dFswing_dwSys(:)];

            obj = obj.writeBlockRhs(r,obj.numLoadCompEqns, ...
                obj.blockJx(obj.map.loadComp,r,Jvals) - F);
            VphA = obj.terminalVphA(obj.loadComp_bus);
            obj = obj.addBlockRhs(r,idxLoadComp_Vs, ...
                reshape(pagemtimes(reshape(dParkV,2,2,n),VphA),2,n));
            obj = obj.addBlockRhs(r,idxLoadComp_WSlip,dFslip_dwSys*ws);
            obj = obj.addBlockRhs(r,idxLoadComp_Is,dFs_dwSys*ws);
            obj = obj.addBlockRhs(r,idxLoadComp_Wr,dFswing_dwSys*ws);
        end
    end

    methods
        %% ================= EMT verification twin ======================
        function emt = buildEmt(obj)
            %BUILDEMT  Construct the time-domain EMT model of this system.
            %
            %   emt = salt.buildEmt() returns an EMT object carrying exactly
            %   the devices and parameters SALT solved, initialised from the
            %   SALT solution, with its stamp maps built and its rotor
            %   angles wrapped. Integrating it forward should leave the
            %   solution where it started - that is the verification.

            emt = EMT(obj.deltaT,obj.tSim(1),obj.tSim(end), ...
                obj.baseMVA,obj.bus);
            emt.bus_BaseKv = obj.bus_BaseKv;
            emt.bus_IBase  = obj.bus_IBase;

            emt = emt.setTxParams(obj.tx_busFrom,obj.tx_busTo, ...
                obj.tx_R,obj.tx_L,obj.tx_C);
            emt = emt.setXfmrParams(obj.xfmr_busFrom,obj.xfmr_busTo, ...
                obj.xfmr_R,obj.xfmr_L);
            if ~isempty(obj.load_bus)
                emt = emt.setLoadParams(obj.load_bus,obj.load_Z);
            end
            emt = emt.setGenrouParams(obj.genrou_bus, ...
                obj.genrou_EsBase,obj.genrou_IsBase,obj.genrou_BaseMVA, ...
                obj.genrou_MW0,obj.genrou_Mvar0, ...
                obj.genrou_Ra,obj.genrou_Rfd,obj.genrou_R1d, ...
                obj.genrou_R1q,obj.genrou_R2q, ...
                obj.genrou_Lad,obj.genrou_Laq,obj.genrou_L0,obj.genrou_Ll, ...
                obj.genrou_Lffd,obj.genrou_Lf1d,obj.genrou_L11d, ...
                obj.genrou_L11q,obj.genrou_L22q, ...
                obj.genrou_H,obj.genrou_KD);
            emt = emt.setExcDc1Params(obj.excDc1_bus,obj.excDc1_Tr, ...
                obj.excDc1_Ta,obj.excDc1_Tc,obj.excDc1_Tb,obj.excDc1_Te, ...
                obj.excDc1_Tf1,obj.excDc1_Kf1,obj.excDc1_Ka,obj.excDc1_Ke);
            emt = emt.setTgov1dParams(obj.tgov1d_bus,obj.tgov1d_T1, ...
                obj.tgov1d_T2,obj.tgov1d_T3,obj.tgov1d_R,obj.tgov1d_Dt);
            if obj.numIbrGFM > 0
                emt = emt.setIbrGFMParams(obj.ibrGFM_bus, ...
                    obj.ibrGFM_BaseKV,obj.ibrGFM_BaseKA,obj.ibrGFM_MVA_base, ...
                    obj.ibrGFM_MW,obj.ibrGFM_MVar, ...
                    obj.ibrGFM_Lf,obj.ibrGFM_Rf,obj.ibrGFM_Cf, ...
                    obj.ibrGFM_Rcap,obj.ibrGFM_Lc,obj.ibrGFM_Rc, ...
                    obj.ibrGFM_kC_i,obj.ibrGFM_kC_p,obj.ibrGFM_GC, ...
                    obj.ibrGFM_kV_i,obj.ibrGFM_kV_p,obj.ibrGFM_GV, ...
                    obj.ibrGFM_wmeas,obj.ibrGFM_droopPercentP, ...
                    obj.ibrGFM_droopPercentQ,obj.ibrGFM_w_VR,obj.ibrGFM_Rv);
            end
            if obj.numLoadComp > 0
                emt = emt.setLoadCompParams(obj.loadComp_bus, ...
                    obj.loadComp_VBase,obj.loadComp_IBase, ...
                    obj.loadComp_MVA_base, ...
                    obj.loadComp_MW,obj.loadComp_MVar, ...
                    obj.loadCompIM_powerPercent, ...
                    obj.loadCompIM_Rs,obj.loadCompIM_Rr, ...
                    obj.loadCompIM_Lss,obj.loadCompIM_Lrr,obj.loadCompIM_Lm, ...
                    obj.loadCompIM_H,obj.loadCompIM_m, ...
                    obj.loadCompZip_PQabc,obj.loadCompZip_Tau);
            end

            % Seed the time-domain state from the SALT phasor solution
            emt = setEmtFromSalt(emt,obj,emt.getSysSize, ...
                obj.tSim,obj.deltaT,obj.x_v);
            emt = emt.buildStampMaps;
            emt.boundThetaRadiansInitializations;
        end
    end

end
