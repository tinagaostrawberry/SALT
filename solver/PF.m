classdef PF < handle
%PF  Positive-sequence power flow in the equivalent-circuit formulation.
%
%   [pfResults,w] = PF.solve(mpc,settings)
%
%   mpc       A MATPOWER case struct. The solver works in EXTERNAL bus
%             numbering from end to end - there is no ext2int step - so
%             every bus number that goes in comes back out unchanged.
%   settings  Struct of solver options:
%
%               Tolerance          Newton relative tolerance.
%               Max_Iters          Newton iteration cap.
%               Limiting           .UseLimiting, .Limit - clamp the size of
%                                  each Newton step, element by element.
%               ZipLoad            .UseZip and, when it is set, the (Z,I,P)
%                                  weight matrices .pw_IM/.qw_IM/.pw_ZIP/
%                                  .qw_ZIP (3 x nBus) against .bus, plus
%                                  .IMpowerPercent/.IMpowerPercent_react,
%                                  the motor share of each load.
%               FreqDeviat         .UseFreqDeviat and, when it is set, the
%                                  droop tables .dP_dw_gen, .dP_dw_load and
%                                  .dQ_dw_load, each [busNum sensitivity].
%
%   pfResults is a copy of mpc with bus VM/VA (pu and DEGREES), bus PD/QD
%   and gen PG/QG written back in MW/MVAr, plus .success, .iterations and
%   .computeTime. w is the converged system frequency in rad/s, or NaN when
%   frequency is not a variable.
%
%   Rather than the equations MATPOWER solves, this is the current-injection
%   ("equivalent circuit") formulation: the unknowns are the real and
%   imaginary parts of each bus voltage, so every device contributes a
%   current-balance row pair and the network itself is linear. PV buses add
%   a reactive-power unknown; the slack adds either two injected-current
%   unknowns or, when frequency is a variable, a Q unknown and one global
%   frequency unknown. See CALCPFNODEINDICES for the layout.
%
%   The Jacobian is never written element by element. Every stamp is emitted
%   as a (value,row,col) triplet - the row and column maps are built once by
%   PREPARESTAMPING, and the values are recomputed per Newton iteration
%   across all devices of a kind at once by the STAMPOFPF*_SS kernels.

    properties (Constant)
        % Absolute floor on the convergence test. A node whose value barely
        % moves in absolute terms is converged even if its relative change
        % is large, which matters for the unknowns that sit near zero.
        convAbsTol = 1e-8
    end

    properties
        %% General
        mpc                           % case the solve was built from
        baseMVA
        w0 = 2*pi*60                  % nominal system frequency (rad/s)

        tol
        max_iters
        enable_limiting
        limit
        useFreqDeviat

        NR_count
        numNodes                      % length of the unknown vector
        v_init

        %% Buses
        % The bus map and the node counter are instance state. The original
        % Buses class kept both in persistent variables, which is why it
        % needed a "clear Buses" before every run and still returned stale
        % indices when that was forgotten.
        numBus = 0
        bus_num
        bus_type
        bus_Vm                        % initial magnitude, pu
        bus_Va                        % initial angle, DEGREES
        bus_node_Vr
        bus_node_Vi
        bus_node_Q                    % 0 where the bus has no Q unknown

        %% Slack
        slack_bus
        slack_Vset                    % [] when frequency is a variable
        slack_AngSet                  % degrees; [] likewise
        slack_imScalConst             % imag/real ratio pinning the angle
        slack_node_Ir
        slack_node_Ii
        node_FreqSlack                % global frequency unknown
        idxW                          % alias of node_FreqSlack, [] if unused

        %% Branches (PI line)
        numBranch = 0
        branch_busFrom
        branch_busTo
        branch_r
        branch_x
        branch_b
        branch_status

        %% Transformers (PI line behind a complex tap)
        numXfmr = 0
        xfmr_busFrom
        xfmr_busTo
        xfmr_r
        xfmr_x
        xfmr_b
        xfmr_status
        xfmr_tap
        xfmr_shift                    % degrees

        %% Shunts
        numShunt = 0
        shunt_bus
        shunt_Gs                      % pu admittance
        shunt_Bs
        shunt_status

        %% Generators (PV buses, plus the slack machine when it is dispatched)
        numGen = 0
        gen_bus
        gen_P                         % pu
        gen_Vset
        gen_Qinit                     % pu, LOAD sign convention
        % Limits are parsed but not enforced: this solver never switches a
        % PV bus to PQ on a Q violation.
        gen_Qmax                      % pu, load sign convention
        gen_Qmin
        gen_status
        gen_dPdw
        gen_Pval                      % dispatched P at the solution, pu

        %% Loads
        numLoad = 0
        load_bus
        load_P                        % pu
        load_Q
        load_status
        load_UseZip
        load_pw_IM                    % [numLoad x 3] (Z,I,P) weights
        load_qw_IM
        load_pw_ZIP
        load_qw_ZIP
        load_IMpowerPercent           % motor share of P
        load_IMpowerPercentReact      % motor share of Q
        load_dPdw
        load_dQdw
        % ZIP coefficients, derived at the seeded voltage (see
        % PFVALINITIALIZATIONS) rather than parsed
        load_Pz_IM,  load_Pi_IM,  load_Pp_IM
        load_Qz_IM,  load_Qi_IM,  load_Qp_IM
        load_Pz_ZIP, load_Pi_ZIP, load_Pp_ZIP
        load_Qz_ZIP, load_Qi_ZIP, load_Qp_ZIP
        load_Pval                     % demand at the solution, pu
        load_Qval

        %% Triplet stamping
        % Row/col maps are built once by PREPARESTAMPING. Their leading
        % numel(stampValsLinear) entries are the linear devices and never
        % change; the rest are the nonlinear devices, recomputed per
        % iteration. The RHS works the same way but accumulates through
        % ACCUMARRAY instead of SPARSE.
        stampIdxRow
        stampJdxRow
        stampValsLinear
        I_idxRow
        I_valsLinear

        % Per-kind gathers of the active devices, so a Newton iteration
        % touches flat arrays only.
        numGenAct = 0
        genAct_idx                    % rows of gen_* these came from
        genAct_idxVr, genAct_idxVi, genAct_idxQ
        genAct_P, genAct_Vset, genAct_dPdw

        numLoadZip = 0
        loadZip_idx
        loadZip_idxVr, loadZip_idxVi
        loadZip_Pz_ZIP, loadZip_Pi_ZIP, loadZip_Pp_ZIP
        loadZip_Qz_ZIP, loadZip_Qi_ZIP, loadZip_Qp_ZIP
        loadZip_Pz_IM,  loadZip_Pi_IM,  loadZip_Pp_IM
        loadZip_Qz_IM,  loadZip_Qi_IM,  loadZip_Qp_IM
        loadZip_dPdw, loadZip_dQdw

        numLoadStd = 0
        loadStd_idx
        loadStd_idxVr, loadStd_idxVi
        loadStd_P, loadStd_Q
        loadStd_IMpowerPercent, loadStd_IMpowerPercentReact
        loadStd_dPdw, loadStd_dQdw
    end

    methods (Static)
        %% ========================= Entry point ========================
        function [pfResults,w] = solve(mpc,settings)
            %SOLVE  Run the power flow and hand back a MATPOWER case copy.
            pf = PF(mpc,settings);
            [v,converge,timeSim] = pf.runPowerFlow();
            [pfResults,w] = pf.processResults(v,converge,timeSim);
        end
    end

    methods
        %% ============================ Class ===========================
        function obj = PF(mpc,settings)
            c = PF.matpowerConstants;

            obj.tol             = settings.Tolerance;
            obj.max_iters       = settings.Max_Iters;
            obj.enable_limiting = settings.Limiting.UseLimiting;
            obj.limit           = settings.Limiting.Limit;
            obj.useFreqDeviat   = settings.FreqDeviat.UseFreqDeviat;

            obj.baseMVA = mpc.baseMVA;

            % ---- BUG FIX (scriptSetupPF.m:178-191) --------------------
            % A PV bus whose generators are all out of service has no way to
            % hold its voltage, so it has to become a PQ bus. The original
            % demoted it only AFTER building the Generators objects, so the
            % offline machine still got an object - and that constructor
            % leaves Bus empty when status is 0, which made assign_nodes
            % fail on an empty bus number. Demoting first means the bus
            % never claims a Q unknown and the dead machine is never built.
            isPvBus = (mpc.bus(:,c.BUS_TYPE) == c.PV);
            busHasOnlineGen = ismember(mpc.bus(:,c.BUS_I), ...
                mpc.gen(mpc.gen(:,c.GEN_STATUS) > 0, c.GEN_BUS));
            mpc.bus(isPvBus & ~busHasOnlineGen, c.BUS_TYPE) = c.PQ;

            % Kept post-demotion, so the returned case reports the bus types
            % the solve actually used.
            obj.mpc = mpc;

            % ---- Devices ---------------------------------------------
            obj = obj.setBusParams(mpc.bus(:,c.BUS_I),mpc.bus(:,c.BUS_TYPE), ...
                mpc.bus(:,c.VM),mpc.bus(:,c.VA));
            obj = obj.setShuntParamsFromCase(mpc,c);
            obj = obj.setSlackParamsFromCase(mpc,c);
            obj = obj.setGenParamsFromCase(mpc,settings,c);
            obj = obj.setNetworkParamsFromCase(mpc,c);
            obj = obj.setLoadParamsFromCase(mpc,settings,c);

            % ---- Node layout, initial guess, triplet maps -------------
            [obj.bus_node_Vr,obj.bus_node_Vi,obj.bus_node_Q, ...
                obj.slack_node_Ir,obj.slack_node_Ii,obj.node_FreqSlack, ...
                obj.numNodes] = ...
                calcPfNodeIndices(obj.bus_type,obj.useFreqDeviat);
            if obj.useFreqDeviat
                obj.idxW = obj.node_FreqSlack;
            else
                obj.idxW = [];
            end

            [obj,v0] = pfValInitializations(obj);
            obj.v_init = v0;

            obj = obj.prepareStamping;
        end

        %% ======================= Set device params ====================

        % ***** BUSES *****
        function obj = setBusParams(obj,bus_num,bus_type,bus_Vm,bus_Va)
            obj.bus_num  = bus_num(:);
            obj.bus_type = bus_type(:);
            obj.bus_Vm   = bus_Vm(:);
            obj.bus_Va   = bus_Va(:);
            obj.numBus   = numel(obj.bus_num);
            assert(numel(unique(obj.bus_num))==obj.numBus, ...
                'Bus numbers must be unique')
        end

        function idx = getBusIdx(obj,busNums)
            %GETBUSIDX  Bus number -> row of the bus arrays.
            % This replaces the persistent containers.Map the old Buses
            % class carried, so the lookup dies with the object.
            b = busNums(:);
            [ok,idx] = ismember(b,obj.bus_num);
            assert(all(ok),'Bus number not present in the case: %s', ...
                mat2str(b(~ok).'))
        end

        % ***** BRANCHES *****
        function obj = setBranchParams(obj,busFrom,busTo,r,x,b,status)
            obj.branch_busFrom = busFrom(:);
            obj.branch_busTo   = busTo(:);
            obj.branch_r       = r(:);
            obj.branch_x       = x(:);
            obj.branch_b       = b(:);
            obj.branch_status  = logical(status(:));
            obj.numBranch      = numel(obj.branch_busFrom);
        end

        % ***** TRANSFORMERS *****
        function obj = setXfmrParams(obj,busFrom,busTo,r,x,b,status,tap,shift)
            obj.xfmr_busFrom = busFrom(:);
            obj.xfmr_busTo   = busTo(:);
            obj.xfmr_r       = r(:);
            obj.xfmr_x       = x(:);
            obj.xfmr_b       = b(:);
            obj.xfmr_status  = logical(status(:));
            % MATPOWER writes a ratio of 0 to mean "no tap", i.e. 1.0
            tap = tap(:);
            tap(tap==0) = 1;
            obj.xfmr_tap   = tap;
            obj.xfmr_shift = shift(:);
            obj.numXfmr    = numel(obj.xfmr_busFrom);
        end

        % ***** SHUNTS *****
        function obj = setShuntParams(obj,bus,Gs_MW,Bs_MVAr,status)
            obj.shunt_bus    = bus(:);
            % MW/MVAr at 1 pu volts is numerically the pu admittance
            obj.shunt_Gs     = Gs_MW(:)   / obj.baseMVA;
            obj.shunt_Bs     = Bs_MVAr(:) / obj.baseMVA;
            obj.shunt_status = logical(status(:));
            obj.numShunt     = numel(obj.shunt_bus);
        end

        % ***** SLACK *****
        function obj = setSlackParams(obj,bus,Vset,AngSet,imScalConst)
            obj.slack_bus         = bus;
            obj.slack_Vset        = Vset;
            obj.slack_AngSet      = AngSet;
            obj.slack_imScalConst = imScalConst;
        end

        % ***** GENERATORS *****
        function obj = setGenParams(obj,bus,P_MW,Vset, ...
                Qmax_MVAr,Qmin_MVAr,Qinit_MVAr,status,dPdw)
            obj.gen_bus  = bus(:);
            obj.gen_P    = P_MW(:) / obj.baseMVA;
            obj.gen_Vset = Vset(:);
            % Q is carried with the LOAD sign convention, because the
            % unknown these machines share with the bus is the reactive
            % power drawn from the network, not injected into it. That is
            % also why the max/min limits swap places.
            obj.gen_Qinit = -Qinit_MVAr(:) / obj.baseMVA;
            obj.gen_Qmax  = -Qmin_MVAr(:)  / obj.baseMVA;
            obj.gen_Qmin  = -Qmax_MVAr(:)  / obj.baseMVA;
            obj.gen_status = logical(status(:));
            obj.gen_dPdw   = dPdw(:);
            obj.numGen     = numel(obj.gen_bus);
            obj.gen_Pval   = obj.gen_P;
        end

        % ***** LOADS *****
        function obj = setLoadParams(obj,bus,P_MW,Q_MVAr,status,useZip, ...
                pw_IM,qw_IM,pw_ZIP,qw_ZIP, ...
                IMpowerPercent,IMpowerPercentReact,dPdw,dQdw)
            obj.load_bus    = bus(:);
            obj.load_P      = P_MW(:)   / obj.baseMVA;
            obj.load_Q      = Q_MVAr(:) / obj.baseMVA;
            obj.load_status = logical(status(:));
            obj.load_UseZip = logical(useZip(:));
            obj.load_pw_IM  = pw_IM;
            obj.load_qw_IM  = qw_IM;
            obj.load_pw_ZIP = pw_ZIP;
            obj.load_qw_ZIP = qw_ZIP;
            obj.load_IMpowerPercent      = IMpowerPercent(:);
            obj.load_IMpowerPercentReact = IMpowerPercentReact(:);
            obj.load_dPdw = dPdw(:);
            obj.load_dQdw = dQdw(:);
            obj.numLoad   = numel(obj.load_bus);
            % Seeded with the nameplate demand so a load that never gets
            % stamped (out of service) still reports something sane.
            obj.load_Pval = obj.load_P;
            obj.load_Qval = obj.load_Q;
        end
    end

    methods (Hidden)
        %% ================= Case data -> device params =================
        function obj = setShuntParamsFromCase(obj,mpc,c)
            % An isolated bus (type 4) keeps its shunt out of the network.
            idx = find((mpc.bus(:,c.GS)~=0 | mpc.bus(:,c.BS)~=0) & ...
                mpc.bus(:,c.BUS_TYPE)~=4);
            obj = obj.setShuntParams(mpc.bus(idx,c.BUS_I), ...
                mpc.bus(idx,c.GS),mpc.bus(idx,c.BS),true(numel(idx),1));
        end

        function obj = setSlackParamsFromCase(obj,mpc,c)
            idxRef = find(mpc.bus(:,c.BUS_TYPE)==c.REF);
            assert(isscalar(idxRef),'PF expects exactly one slack bus')
            busNum = mpc.bus(idxRef,c.BUS_I);

            if obj.useFreqDeviat
                % With frequency free, the slack no longer fixes its
                % voltage - it only pins the angle - so the magnitude set
                % point is meaningless and the imag/real ratio takes over.
                % The angle is held by forcing Vi/Vr to the constant ratio
                % of the seed voltage, since angle(Vr+1j*Vi) is unchanged
                % when both parts are scaled by the same factor.
                % (Inlined from the old staticDevices_getSlackParams.)
                Vslack = mpc.bus(idxRef,c.VM) * ...
                    exp(1j*deg2rad(mpc.bus(idxRef,c.VA)));
                assert(real(Vslack)~=0, ...
                    'PF: slack voltage has zero real part, angle cannot be pinned')
                obj = obj.setSlackParams(busNum,[],[], ...
                    imag(Vslack)/real(Vslack));
            else
                obj = obj.setSlackParams(busNum, ...
                    mpc.bus(idxRef,c.VM),mpc.bus(idxRef,c.VA),[]);
            end
        end

        function obj = setGenParamsFromCase(obj,mpc,settings,c)
            assert(numel(unique(mpc.gen(:,c.GEN_BUS)))==size(mpc.gen,1), ...
                'PF expects at most one generator per bus')

            genRows = zeros(0,1);
            if obj.useFreqDeviat
                % The slack machine is dispatched like any other generator
                % once frequency is free, so it needs a device of its own.
                [ok,iSlack] = ismember(obj.slack_bus,mpc.gen(:,c.GEN_BUS));
                assert(ok,'UseFreqDeviat needs a generator at the slack bus')
                genRows = iSlack;
            end
            isPvGen = ismember(mpc.gen(:,c.GEN_BUS), ...
                mpc.bus(mpc.bus(:,c.BUS_TYPE)==c.PV,c.BUS_I));
            genRows = [genRows; find(isPvGen)];

            % Out-of-service machines stamp nothing and report nothing (the
            % result write-back zeroes their PG/QG regardless), so they are
            % dropped here instead of being built as objects that could not
            % be given a bus. See the demotion fix in the constructor.
            genRows = genRows(mpc.gen(genRows,c.GEN_STATUS) > 0);

            if obj.useFreqDeviat
                tbl = settings.FreqDeviat.dP_dw_gen;
                [ok,loc] = ismember(mpc.gen(genRows,c.GEN_BUS),tbl(:,1));
                assert(all(ok), ...
                    'FreqDeviat.dP_dw_gen has no row for generator bus(es): %s', ...
                    mat2str(mpc.gen(genRows(~ok),c.GEN_BUS).'))
                dPdw = tbl(loc,2);
            else
                dPdw = zeros(numel(genRows),1);
            end

            obj = obj.setGenParams(mpc.gen(genRows,c.GEN_BUS), ...
                mpc.gen(genRows,c.PG),mpc.gen(genRows,c.VG), ...
                mpc.gen(genRows,c.QMAX),mpc.gen(genRows,c.QMIN), ...
                mpc.gen(genRows,c.QG), ...
                mpc.gen(genRows,c.GEN_STATUS),dPdw);
        end

        function obj = setNetworkParamsFromCase(obj,mpc,c)
            % A tap that is neither 0 nor 1, or any phase shift, makes the
            % branch a transformer; everything else is a plain PI line.
            isXfmr = (mpc.branch(:,c.TAP)~=0 & mpc.branch(:,c.TAP)~=1) | ...
                (mpc.branch(:,c.SHIFT)~=0);

            iBr = find(~isXfmr);
            obj = obj.setBranchParams(mpc.branch(iBr,c.F_BUS), ...
                mpc.branch(iBr,c.T_BUS),mpc.branch(iBr,c.BR_R), ...
                mpc.branch(iBr,c.BR_X),mpc.branch(iBr,c.BR_B), ...
                mpc.branch(iBr,c.BR_STATUS));

            iXf = find(isXfmr);
            obj = obj.setXfmrParams(mpc.branch(iXf,c.F_BUS), ...
                mpc.branch(iXf,c.T_BUS),mpc.branch(iXf,c.BR_R), ...
                mpc.branch(iXf,c.BR_X),mpc.branch(iXf,c.BR_B), ...
                mpc.branch(iXf,c.BR_STATUS), ...
                mpc.branch(iXf,c.TAP),mpc.branch(iXf,c.SHIFT));
        end

        function obj = setLoadParamsFromCase(obj,mpc,settings,c)
            idx = find(mpc.bus(:,c.PD)~=0 | mpc.bus(:,c.QD)~=0);
            n = numel(idx);

            busL   = mpc.bus(idx,c.BUS_I);
            P      = mpc.bus(idx,c.PD);
            Q      = mpc.bus(idx,c.QD);
            status = (mpc.bus(idx,c.BUS_TYPE) ~= 4);

            pw_IM = zeros(n,3); qw_IM = zeros(n,3);
            pw_ZIP = zeros(n,3); qw_ZIP = zeros(n,3);
            IMp = zeros(n,1); IMpReact = zeros(n,1);
            useZip = false(n,1);
            dPdw = zeros(n,1); dQdw = zeros(n,1);

            if n > 0
                zl = settings.ZipLoad;
                useZip(:) = logical(zl.UseZip);

                % ZipLoad.bus indexes the per-bus ZIP columns. A load bus
                % missing from it simply has no motor share.
                % (ZipLoad.pw / ZipLoad.qw are deliberately never read - the
                % solver only uses the IM/ZIP split.)
                [inZip,locZip] = ismember(busL,zl.bus(:));
                IMp(inZip)      = zl.IMpowerPercent(locZip(inZip));
                IMpReact(inZip) = zl.IMpowerPercent_react(locZip(inZip));

                if zl.UseZip
                    assert(all(inZip), ...
                        'ZipLoad.bus must cover every load bus when UseZip is true: %s', ...
                        mat2str(busL(~inZip).'))
                    pw_IM  = zl.pw_IM(:,locZip).';
                    qw_IM  = zl.qw_IM(:,locZip).';
                    pw_ZIP = zl.pw_ZIP(:,locZip).';
                    qw_ZIP = zl.qw_ZIP(:,locZip).';
                    % Motor and static shares must together account for the
                    % whole load, or the ZIP split silently rescales it.
                    assert(all(abs(sum([pw_IM pw_ZIP],2)-1) < 1e-6))
                    assert(all(abs(sum([qw_IM qw_ZIP],2)-1) < 1e-6))
                end

                if obj.useFreqDeviat
                    [okP,locP] = ismember(busL, ...
                        settings.FreqDeviat.dP_dw_load(:,1));
                    assert(all(okP), ...
                        'FreqDeviat.dP_dw_load has no row for load bus(es): %s', ...
                        mat2str(busL(~okP).'))
                    dPdw = settings.FreqDeviat.dP_dw_load(locP,2);

                    [okQ,locQ] = ismember(busL, ...
                        settings.FreqDeviat.dQ_dw_load(:,1));
                    assert(all(okQ), ...
                        'FreqDeviat.dQ_dw_load has no row for load bus(es): %s', ...
                        mat2str(busL(~okQ).'))
                    dQdw = settings.FreqDeviat.dQ_dw_load(locQ,2);
                end
            end

            obj = obj.setLoadParams(busL,P,Q,status,useZip, ...
                pw_IM,qw_IM,pw_ZIP,qw_ZIP,IMp,IMpReact,dPdw,dQdw);
        end
    end

    methods
        %% ======================= Triplet stamping =====================
        function obj = prepareStamping(obj)
            %PREPARESTAMPING  Build the fixed row/col maps and the constant
            %   linear values, and gather the active devices into the flat
            %   arrays the per-iteration kernels consume.

            % ---- Linear devices, stamped once -------------------------
            brOn = obj.branch_status;
            xfOn = obj.xfmr_status;
            shOn = obj.shunt_status;

            brN = obj.getBusIdx(obj.branch_busFrom(brOn));
            brM = obj.getBusIdx(obj.branch_busTo(brOn));
            xfN = obj.getBusIdx(obj.xfmr_busFrom(xfOn));
            xfM = obj.getBusIdx(obj.xfmr_busTo(xfOn));
            shI = obj.getBusIdx(obj.shunt_bus(shOn));
            slI = obj.getBusIdx(obj.slack_bus);

            [yRow,yCol,yVal,iRow,iVal] = stampOfPfLinearMatrices_SS( ...
                obj.bus_node_Vr(brN),obj.bus_node_Vi(brN), ...
                obj.bus_node_Vr(brM),obj.bus_node_Vi(brM), ...
                obj.branch_r(brOn),obj.branch_x(brOn),obj.branch_b(brOn), ...
                ...
                obj.bus_node_Vr(xfN),obj.bus_node_Vi(xfN), ...
                obj.bus_node_Vr(xfM),obj.bus_node_Vi(xfM), ...
                obj.xfmr_r(xfOn),obj.xfmr_x(xfOn),obj.xfmr_b(xfOn), ...
                obj.xfmr_tap(xfOn),obj.xfmr_shift(xfOn), ...
                ...
                obj.bus_node_Vr(shI),obj.bus_node_Vi(shI), ...
                obj.shunt_Gs(shOn),obj.shunt_Bs(shOn), ...
                ...
                obj.bus_node_Vr(slI),obj.bus_node_Vi(slI), ...
                obj.slack_node_Ir,obj.slack_node_Ii, ...
                obj.slack_Vset,obj.slack_AngSet, ...
                obj.slack_imScalConst,obj.node_FreqSlack, ...
                ...
                obj.useFreqDeviat);

            obj.stampValsLinear = yVal;
            obj.I_valsLinear    = iVal;
            idxRow   = yRow;
            idxCol   = yCol;
            idxRow_I = iRow;

            % ---- Gather the active generators -------------------------
            act = obj.gen_status;
            obj.genAct_idx = find(act);
            ibG = obj.getBusIdx(obj.gen_bus(act));
            obj.genAct_idxVr = obj.bus_node_Vr(ibG);
            obj.genAct_idxVi = obj.bus_node_Vi(ibG);
            obj.genAct_idxQ  = obj.bus_node_Q(ibG);
            obj.genAct_P     = obj.gen_P(act);
            obj.genAct_Vset  = obj.gen_Vset(act);
            obj.genAct_dPdw  = obj.gen_dPdw(act);
            obj.numGenAct    = numel(obj.genAct_idx);
            % Every generator must sit on a bus that owns a Q unknown; this
            % is what the old Generators.assign_nodes check was guarding.
            assert(all(obj.genAct_idxQ > 0), ...
                'Generator on a bus with no reactive-power unknown')

            % ---- Gather the active loads, split by load model ---------
            actL    = obj.load_status;
            zipMask = actL &  obj.load_UseZip;
            stdMask = actL & ~obj.load_UseZip;

            obj.loadZip_idx = find(zipMask);
            ibZ = obj.getBusIdx(obj.load_bus(zipMask));
            obj.loadZip_idxVr = obj.bus_node_Vr(ibZ);
            obj.loadZip_idxVi = obj.bus_node_Vi(ibZ);
            obj.loadZip_Pz_ZIP = obj.load_Pz_ZIP(zipMask);
            obj.loadZip_Pi_ZIP = obj.load_Pi_ZIP(zipMask);
            obj.loadZip_Pp_ZIP = obj.load_Pp_ZIP(zipMask);
            obj.loadZip_Qz_ZIP = obj.load_Qz_ZIP(zipMask);
            obj.loadZip_Qi_ZIP = obj.load_Qi_ZIP(zipMask);
            obj.loadZip_Qp_ZIP = obj.load_Qp_ZIP(zipMask);
            obj.loadZip_Pz_IM  = obj.load_Pz_IM(zipMask);
            obj.loadZip_Pi_IM  = obj.load_Pi_IM(zipMask);
            obj.loadZip_Pp_IM  = obj.load_Pp_IM(zipMask);
            obj.loadZip_Qz_IM  = obj.load_Qz_IM(zipMask);
            obj.loadZip_Qi_IM  = obj.load_Qi_IM(zipMask);
            obj.loadZip_Qp_IM  = obj.load_Qp_IM(zipMask);
            obj.loadZip_dPdw   = obj.load_dPdw(zipMask);
            obj.loadZip_dQdw   = obj.load_dQdw(zipMask);
            obj.numLoadZip     = numel(obj.loadZip_idx);

            obj.loadStd_idx = find(stdMask);
            ibS = obj.getBusIdx(obj.load_bus(stdMask));
            obj.loadStd_idxVr = obj.bus_node_Vr(ibS);
            obj.loadStd_idxVi = obj.bus_node_Vi(ibS);
            obj.loadStd_P = obj.load_P(stdMask);
            obj.loadStd_Q = obj.load_Q(stdMask);
            obj.loadStd_IMpowerPercent      = obj.load_IMpowerPercent(stdMask);
            obj.loadStd_IMpowerPercentReact = obj.load_IMpowerPercentReact(stdMask);
            obj.loadStd_dPdw = obj.load_dPdw(stdMask);
            obj.loadStd_dQdw = obj.load_dQdw(stdMask);
            obj.numLoadStd   = numel(obj.loadStd_idx);

            % ---- Nonlinear sparsity pattern (row/col only) ------------
            % Order MUST match the value order the STAMPOFPF*_SS kernels
            % return, since the two are zipped together positionally.
            nG = obj.numGenAct;
            if nG > 0
                gVr = obj.genAct_idxVr; gVi = obj.genAct_idxVi;
                gQ  = obj.genAct_idxQ;
                rowsG = [gVr; gVr; gVr; gVi; gVi; gVi; gQ;  gQ];
                colsG = [gQ;  gVr; gVi; gQ;  gVr; gVi; gVr; gVi];
                if obj.useFreqDeviat
                    rowsG = [rowsG; gVr; gVi];
                    colsG = [colsG; repmat(obj.idxW,nG,1); repmat(obj.idxW,nG,1)];
                end
                idxRow   = [idxRow; rowsG];
                idxCol   = [idxCol; colsG];
                idxRow_I = [idxRow_I; gVr; gVi; gQ];
            end

            nZ = obj.numLoadZip;
            if nZ > 0
                zVr = obj.loadZip_idxVr; zVi = obj.loadZip_idxVi;
                rowsZ = [zVr; zVr; zVi; zVi];
                colsZ = [zVr; zVi; zVr; zVi];
                if obj.useFreqDeviat
                    rowsZ = [rowsZ; zVr; zVi];
                    colsZ = [colsZ; repmat(obj.idxW,nZ,1); repmat(obj.idxW,nZ,1)];
                end
                idxRow   = [idxRow; rowsZ];
                idxCol   = [idxCol; colsZ];
                idxRow_I = [idxRow_I; zVr; zVi];
            end

            nS = obj.numLoadStd;
            if nS > 0
                sVr = obj.loadStd_idxVr; sVi = obj.loadStd_idxVi;
                rowsS = [sVr; sVr; sVi; sVi];
                colsS = [sVr; sVi; sVr; sVi];
                if obj.useFreqDeviat
                    rowsS = [rowsS; sVr; sVi];
                    colsS = [colsS; repmat(obj.idxW,nS,1); repmat(obj.idxW,nS,1)];
                end
                idxRow   = [idxRow; rowsS];
                idxCol   = [idxCol; colsS];
                idxRow_I = [idxRow_I; sVr; sVi];
            end

            obj.stampIdxRow = idxRow;
            obj.stampJdxRow = idxCol;
            obj.I_idxRow    = idxRow_I;
        end

        function [stampVals,iVals] = assembleStampValues(obj,v,stampVals,iVals)
            %ASSEMBLESTAMPVALUES  Refresh the nonlinear part of the value
            %   arrays at the current guess. The linear head of both arrays
            %   is already in place and is left untouched.
            offsetY = numel(obj.stampValsLinear);
            offsetI = numel(obj.I_valsLinear);

            if obj.useFreqDeviat
                w = v(obj.idxW);
            else
                w = [];
            end

            if obj.numGenAct > 0
                [yV,iV] = stampOfPfGenMatrices_SS( ...
                    v(obj.genAct_idxVr),v(obj.genAct_idxVi), ...
                    v(obj.genAct_idxQ),w, ...
                    obj.genAct_P,obj.genAct_Vset,obj.genAct_dPdw, ...
                    obj.w0,obj.useFreqDeviat);
                stampVals(offsetY+(1:numel(yV))) = yV;
                offsetY = offsetY + numel(yV);
                iVals(offsetI+(1:numel(iV))) = iV;
                offsetI = offsetI + numel(iV);
            end

            if obj.numLoadZip > 0
                [yV,iV] = stampOfPfLoadZipMatrices_SS( ...
                    v(obj.loadZip_idxVr),v(obj.loadZip_idxVi),w, ...
                    obj.loadZip_Pz_ZIP,obj.loadZip_Pi_ZIP,obj.loadZip_Pp_ZIP, ...
                    obj.loadZip_Qz_ZIP,obj.loadZip_Qi_ZIP,obj.loadZip_Qp_ZIP, ...
                    obj.loadZip_Pz_IM,obj.loadZip_Pi_IM,obj.loadZip_Pp_IM, ...
                    obj.loadZip_Qz_IM,obj.loadZip_Qi_IM,obj.loadZip_Qp_IM, ...
                    obj.loadZip_dPdw,obj.loadZip_dQdw, ...
                    obj.w0,obj.useFreqDeviat);
                stampVals(offsetY+(1:numel(yV))) = yV;
                offsetY = offsetY + numel(yV);
                iVals(offsetI+(1:numel(iV))) = iV;
                offsetI = offsetI + numel(iV);
            end

            if obj.numLoadStd > 0
                [yV,iV] = stampOfPfLoadStdMatrices_SS( ...
                    v(obj.loadStd_idxVr),v(obj.loadStd_idxVi),w, ...
                    obj.loadStd_P,obj.loadStd_Q, ...
                    obj.loadStd_IMpowerPercent, ...
                    obj.loadStd_IMpowerPercentReact, ...
                    obj.loadStd_dPdw,obj.loadStd_dQdw, ...
                    obj.w0,obj.useFreqDeviat);
                stampVals(offsetY+(1:numel(yV))) = yV;
                iVals(offsetI+(1:numel(iV))) = iV;
            end
        end

        %% =========================== Solve ============================
        function [v,converge,timeSim] = runPowerFlow(obj)
            %RUNPOWERFLOW  Newton-Raphson on the current-injection system.

            v = obj.v_init;
            n = obj.numNodes;

            % The linear head of both value arrays is constant, so it is
            % written once and only the nonlinear tail is refreshed below.
            stampVals = zeros(numel(obj.stampIdxRow),1);
            stampVals(1:numel(obj.stampValsLinear)) = obj.stampValsLinear;
            iVals = zeros(numel(obj.I_idxRow),1);
            iVals(1:numel(obj.I_valsLinear)) = obj.I_valsLinear;

            obj.NR_count = 0;
            converge = false;

            TIMER = tic;
            while (obj.NR_count < obj.max_iters) && ~converge
                [stampVals,iVals] = obj.assembleStampValues(v,stampVals,iVals);

                Y = sparse(obj.stampIdxRow,obj.stampJdxRow,stampVals,n,n);
                I = accumarray(obj.I_idxRow,iVals,[n,1]);

                v_prev = v;

                % Solves Y*v = I. The original called KLU here
                % (klu(Y,'\',I)); KLU is not installed, and MATLAB's sparse
                % backslash performs the same solve.
                v = Y \ I;

                % Convergence is judged on the raw Newton step, before any
                % limiting, so a clamped step never reads as converged.
                converge = obj.checkError(v,v_prev);

                if obj.enable_limiting
                    dv = v - v_prev;
                    idxViolat = abs(dv) > obj.limit;
                    dv(idxViolat) = sign(dv(idxViolat))*obj.limit;
                    v = v_prev + dv;
                end

                obj.NR_count = obj.NR_count + 1;
            end
            timeSim = toc(TIMER);

            % Per-device derived state at the CONVERGED point - one Newton
            % step past the last assembly, so it has to be recomputed.
            obj.writeBackDeviceState(v);
        end

        function convergence = checkError(obj,v,v_prev)
            % A node passes on either measure: absolute for the unknowns
            % that sit near zero, relative for everything else.
            absErr = abs(v - v_prev);
            relErr = abs((v - v_prev) ./ (v_prev + eps));
            convergence = all(absErr < obj.convAbsTol | relErr < obj.tol);
        end

        function obj = writeBackDeviceState(obj,v)
            %WRITEBACKDEVICESTATE  Dispatched P (and load P/Q) at the
            %   solution, which differ from the nameplate values whenever
            %   the load is voltage dependent or frequency has moved.
            if obj.useFreqDeviat
                w = v(obj.idxW);
            else
                w = [];
            end

            if obj.numGenAct > 0
                [~,~,P_real] = stampOfPfGenMatrices_SS( ...
                    v(obj.genAct_idxVr),v(obj.genAct_idxVi), ...
                    v(obj.genAct_idxQ),w, ...
                    obj.genAct_P,obj.genAct_Vset,obj.genAct_dPdw, ...
                    obj.w0,obj.useFreqDeviat);
                obj.gen_Pval(obj.genAct_idx) = P_real;
            end

            if obj.numLoadZip > 0
                [~,~,Pval,Qval] = stampOfPfLoadZipMatrices_SS( ...
                    v(obj.loadZip_idxVr),v(obj.loadZip_idxVi),w, ...
                    obj.loadZip_Pz_ZIP,obj.loadZip_Pi_ZIP,obj.loadZip_Pp_ZIP, ...
                    obj.loadZip_Qz_ZIP,obj.loadZip_Qi_ZIP,obj.loadZip_Qp_ZIP, ...
                    obj.loadZip_Pz_IM,obj.loadZip_Pi_IM,obj.loadZip_Pp_IM, ...
                    obj.loadZip_Qz_IM,obj.loadZip_Qi_IM,obj.loadZip_Qp_IM, ...
                    obj.loadZip_dPdw,obj.loadZip_dQdw, ...
                    obj.w0,obj.useFreqDeviat);
                obj.load_Pval(obj.loadZip_idx) = Pval;
                obj.load_Qval(obj.loadZip_idx) = Qval;
            end

            if obj.numLoadStd > 0
                [~,~,Pval,Qval] = stampOfPfLoadStdMatrices_SS( ...
                    v(obj.loadStd_idxVr),v(obj.loadStd_idxVi),w, ...
                    obj.loadStd_P,obj.loadStd_Q, ...
                    obj.loadStd_IMpowerPercent, ...
                    obj.loadStd_IMpowerPercentReact, ...
                    obj.loadStd_dPdw,obj.loadStd_dQdw, ...
                    obj.w0,obj.useFreqDeviat);
                obj.load_Pval(obj.loadStd_idx) = Pval;
                obj.load_Qval(obj.loadStd_idx) = Qval;
            end
        end

        %% ========================== Results ===========================
        function [pfResults,w] = processResults(obj,v,converge,timeSim)
            %PROCESSRESULTS  Write the solution back into a case copy.
            c = PF.matpowerConstants;

            % ---- Bus voltages ----------------------------------------
            Vr = v(obj.bus_node_Vr);
            Vi = v(obj.bus_node_Vi);
            Vm     = sqrt(Vr.^2 + Vi.^2);
            Va_rad = atan2(Vi,Vr);

            % ---- Generator dispatch ----------------------------------
            % Q comes straight off the bus unknown, negated back out of the
            % load sign convention the solver carries it in.
            idxBusGen = obj.getBusIdx(obj.gen_bus);
            Pgen = obj.baseMVA * obj.gen_Pval;
            Qgen = -obj.baseMVA * v(obj.bus_node_Q(idxBusGen));

            % ---- Load demand -----------------------------------------
            % Only voltage- or frequency-dependent loads moved off their
            % nameplate value. The AND on the two droop terms is the
            % original condition, kept as-is.
            useSaved = obj.load_UseZip | ...
                ((obj.load_dPdw~=0) & (obj.load_dQdw~=0));
            Pload = obj.baseMVA * obj.load_P;
            Qload = obj.baseMVA * obj.load_Q;
            Pload(useSaved) = obj.baseMVA * obj.load_Pval(useSaved);
            Qload(useSaved) = obj.baseMVA * obj.load_Qval(useSaved);

            % ---- Assemble the returned case --------------------------
            pfResults = obj.mpc;

            [~,idxBus] = ismember(obj.bus_num,pfResults.bus(:,c.BUS_I));
            pfResults.bus(idxBus,[c.VM c.VA]) = [Vm rad2deg(Va_rad)];

            % BUG FIX (solve.m:29-30): the original also wrote
            %   mpc.gen(:,[VM VA]) = ...
            % but VM/VA are bus columns 8-9, which in the GEN table are
            % GEN_STATUS and PMAX - and it wrote them into a local copy that
            % was never returned. Corrupting-but-discarded; deleted outright.

            [~,idxLoad] = ismember(obj.load_bus,pfResults.bus(:,c.BUS_I));
            pfResults.bus(idxLoad,c.PD) = Pload;
            pfResults.bus(idxLoad,c.QD) = Qload;

            [~,idxGen] = ismember(obj.gen_bus,pfResults.gen(:,c.GEN_BUS));
            pfResults.gen(idxGen,c.PG) = Pgen;
            pfResults.gen(idxGen,c.QG) = Qgen;

            % An out-of-service machine produces nothing, whatever the
            % case file happened to carry.
            idxGenOff = (pfResults.gen(:,c.GEN_STATUS)==0);
            pfResults.gen(idxGenOff,[c.PG c.QG]) = 0;

            % ---- Slack ------------------------------------------------
            % With frequency free the slack machine is an ordinary
            % generator and was already written above; otherwise its output
            % is whatever the injected-current unknowns came out to.
            if ~obj.useFreqDeviat
                idxSlackBus = obj.getBusIdx(obj.slack_bus);
                iComplex = v(obj.slack_node_Ir) + 1j*v(obj.slack_node_Ii);
                vComplex = v(obj.bus_node_Vr(idxSlackBus)) + ...
                    1j*v(obj.bus_node_Vi(idxSlackBus));
                Scomplex = vComplex*conj(iComplex);

                [okSlackGen,idxGenSlack] = ismember(obj.slack_bus, ...
                    pfResults.gen(:,c.GEN_BUS));
                assert(okSlackGen, ...
                    'No generator row at slack bus %g to receive its dispatch', ...
                    obj.slack_bus)
                pfResults.gen(idxGenSlack,[c.PG c.QG]) = ...
                    [real(Scomplex) imag(Scomplex)]*obj.baseMVA;
            end

            % ---- Frequency and solve metadata -------------------------
            if obj.useFreqDeviat
                w = v(obj.node_FreqSlack);
            else
                w = NaN;
            end

            pfResults.success     = converge;
            pfResults.iterations  = obj.NR_count;
            pfResults.computeTime = timeSim;
        end
    end

    methods (Static, Hidden)
        function c = matpowerConstants
            %MATPOWERCONSTANTS  Column indices of the MATPOWER case tables.
            % Taken from MATPOWER itself rather than hard-coded, so a change
            % in the case format is picked up here.
            c = struct;
            [c.PQ,c.PV,c.REF,~,c.BUS_I,c.BUS_TYPE,c.PD,c.QD,c.GS,c.BS, ...
                ~,c.VM,c.VA] = idx_bus;
            [c.GEN_BUS,c.PG,c.QG,c.QMAX,c.QMIN,c.VG,~,c.GEN_STATUS] = idx_gen;
            [c.F_BUS,c.T_BUS,c.BR_R,c.BR_X,c.BR_B,~,~,~, ...
                c.TAP,c.SHIFT,c.BR_STATUS] = idx_brch;
        end
    end
end
