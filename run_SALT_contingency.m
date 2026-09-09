%% Contingency analysis: SALT vs enhanced power flow vs power flow
%
% Sweeps every single-generator outage on the 39-bus system with high
% inverter penetration and counts the limit violations each model predicts.
% Three models see the same contingency:
%
%   SALT     the steady-state EMT solver in solver/SALT.m. Physics-based
%            device models, and the system frequency is an unknown, so a
%            power imbalance moves it off 60 Hz.
%   Enh-PF   PF (solver/PF.m) with the composite loads reduced to per-bus
%            ZIP coefficients and every device given a frequency droop, so
%            the power flow also solves for frequency. Those two data sets
%            are derived FROM the SALT model, by compLoad2Zip_0 and
%            emtDevice2FreqDevInPf.
%   PF       PF with the nominal settings: constant-power loads, fixed
%            60 Hz. This is the conventional answer.
%
% The point of the comparison is which contingencies each model flags:
% bus-voltage limits, generator Q limits and branch ratings.

clear; close all; clc
root = addpathSALT;

%% Case
% The 39-bus network with grid-forming inverters at buses 31-39 and a
% single machine at bus 30. Every device parameter - inverter filter and
% control gains, induction-motor and ZIP data - comes from the workbook.
rawFile    = fullfile(root,'data','case39_modified.m');
deviceFile = fullfile(root,'data','39bus_highPenetrationIBR_SALT.xlsx');
caseName   = "39bus_highPenetrationIBR";

% No load perturbation: the only thing that moves between runs is the
% contingency.
loadScale = 1;

%% Solver settings

% Newton-Raphson, for SALT
totIters_NR  = 50;
relErrTol_NR = 1e-4;
absErrTol_NR = 1e-8;

% Absolute step clamp. Voltages are per unit and the frequency should not
% move far, so no single Newton step should exceed 1.
include_clampLimiting = true;
clampLimit = 1;

% Time step for the sample grid SALT carries for its EMT twin. Nothing here
% simulates in the time domain, but the constructor needs it.
maxDeltaT = 1e-5;

% Power flow. The tolerances are MATPOWER's defaults so the comparison is
% against a conventionally-tuned power flow.
settings_nominal = struct( ...
    'Tolerance',1e-8, ...
    'Max_Iters',10, ...
    'Limiting',struct('UseLimiting',false,'Limit',NaN), ...
    'ZipLoad',struct( ...
        'UseZip',false, ...
        'pw_IM',[],'qw_IM',[],'pw_ZIP',[],'qw_ZIP',[], ...
        'pw',[],'qw',[],'bus',[], ...
        'IMpowerPercent',[],'IMpowerPercent_react',[]), ...
    'FreqDeviat',struct( ...
        'UseFreqDeviat',false, ...
        'dP_dw_gen',[],'dP_dw_load',[],'dQ_dw_load',[]));

%% Base case: build and solve SALT

salt = SALT(rawFile,deviceFile, ...
    maxDeltaT   = maxDeltaT, ...
    loadScale   = loadScale, ...
    lineOutage  = false, ...
    genOutage   = false, ...
    loadOutage  = false, ...
    rerunPfInit = true);
salt.caseName = char(caseName);

% The solved base case. SALT re-runs the power flow to initialise, so
% salt.mpc is already a converged operating point and both solvers start
% from the same one.
scriptSimulate_SALT
mpc_0 = salt.mpc;

bus_VreIm_SALT_0 = salt.x_v(salt.getRowIdxAbc_fromBusNums(salt.bus));
bus_V_SALT_0  = bus_VreIm_SALT_0(1:2:(2*salt.numBus)) + ...
              1j*bus_VreIm_SALT_0(1+(1:2:(2*salt.numBus)));
bus_Vm_SALT_0 = abs(bus_V_SALT_0);
bus_Va_SALT_0 = angle(bus_V_SALT_0);
[I_gen_SALT_0,I_load_SALT_0,I_series_SALT_0] = calcSaltCurrents(salt);

%% Enhanced power-flow settings, derived from the SALT model
% Both conversions read the CONVERGED SALT object, so the enhanced power
% flow is given exactly the load and droop behaviour SALT's device models
% produce at this operating point - not an independently guessed fit.

[Pzip_zip_PF,Qzip_zip_PF,Pzip_IM_PF,Qzip_IM_PF,Pzip_PF,Qzip_PF, ...
    compLoadZipBus_PF,IMpowerPercent_react] = compLoad2Zip(salt);

[dP_dw_gen,dP_dw_ibrGFM,dP_dw_load,dQ_dw_load] = ...
    emtDevice2FreqDevInPf(salt);

settings = settings_nominal;
% ZIP loads, one coefficient triple per bus
settings.ZipLoad.UseZip  = true;
settings.ZipLoad.pw_IM   = Pzip_IM_PF;
settings.ZipLoad.qw_IM   = Qzip_IM_PF;
settings.ZipLoad.pw_ZIP  = Pzip_zip_PF;
settings.ZipLoad.qw_ZIP  = Qzip_zip_PF;
settings.ZipLoad.pw      = Pzip_PF;
settings.ZipLoad.qw      = Qzip_PF;
settings.ZipLoad.bus     = compLoadZipBus_PF;
settings.ZipLoad.IMpowerPercent = salt.loadCompIM_powerPercent(:).';
settings.ZipLoad.IMpowerPercent_react = IMpowerPercent_react;
% Frequency droop. Machines and inverters share one list.
settings.FreqDeviat.UseFreqDeviat = true;
settings.FreqDeviat.dP_dw_gen  = [dP_dw_gen;dP_dw_ibrGFM];
settings.FreqDeviat.dP_dw_load = dP_dw_load;
settings.FreqDeviat.dQ_dw_load = dQ_dw_load;

%% Base case: check all three models agree
% With no contingency the system is at its nominal operating point, so the
% plain power flow, the enhanced power flow and SALT must all return it.
% Any disagreement here means the models were not set up from the same
% data, and every contingency result afterwards would be meaningless.

% MATPOWER column indices. NOTE: do NOT use define_constants here - it
% defines PF as a branch column index (real power at the "from" end), which
% would shadow the PF class for the rest of the script. The idx_* accessors
% give the same numbers without that name.
[~,~,~,~,BUS_I,~,~,~,~,~,~,VM,VA,~,~,VMAX,VMIN] = idx_bus;
[GEN_BUS,~,~,QMAX,QMIN,~,MBASE,~,PMAX,PMIN] = idx_gen;
[F_BUS,T_BUS,~,~,~,RATE_A] = idx_brch;

for idxModel = 1:2
    if idxModel==1
        [mpc_0,w_PF] = PF.solve(mpc_0,settings_nominal);
        % Frequency is not a variable in the nominal settings, so PF
        % reports NaN; the operating point is nominal by construction.
        assert(isnan(w_PF))
        w_PF = salt.w0;
    else
        [mpc_0,w_PF] = PF.solve(mpc_0,settings);
    end
    bus_Vm_pf_0 = mpc_0.bus(:,VM);
    bus_Va_pf_0 = deg2rad(mpc_0.bus(:,VA));
    [I_gen_pf_0,I_load_pf_0,I_series_pf_0] = calcPfCurrents(mpc_0);

    assert(all(abs(bus_Vm_pf_0-bus_Vm_SALT_0)<1e-6))
    assert(all(abs(bus_Va_pf_0-bus_Va_SALT_0)<1e-6))
    assert(all(abs(I_gen_pf_0-I_gen_SALT_0)<1e-6,'all'))
    assert(all(abs(I_load_pf_0-I_load_SALT_0)<1e-6,'all'))
    assert(all(abs(I_series_pf_0-I_series_SALT_0)<1e-6,'all'))
    assert(abs(w_PF-salt.w0)<1e-3)
end

%% Limits to test each contingency against
% Generator limits are on the MACHINE base in the case file and are
% converted to the system (bus) base, because that is the base the currents
% and voltages above are in.

bus     = mpc_0.bus(:,BUS_I);
busVmin = mpc_0.bus(:,[BUS_I VMIN]);
busVmax = mpc_0.bus(:,[BUS_I VMAX]);

genQmax = mpc_0.gen(:,[GEN_BUS QMAX]);
genQmin = mpc_0.gen(:,[GEN_BUS QMIN]);
genPmax = mpc_0.gen(:,[GEN_BUS PMAX]);
genPmin = mpc_0.gen(:,[GEN_BUS PMIN]);
genQmax(:,2) = genQmax(:,2)./mpc_0.gen(:,MBASE);
genQmin(:,2) = genQmin(:,2)./mpc_0.gen(:,MBASE);
genPmax(:,2) = genPmax(:,2)./mpc_0.gen(:,MBASE);
genPmin(:,2) = genPmin(:,2)./mpc_0.gen(:,MBASE);

brchRateMva = mpc_0.branch(:,[F_BUS T_BUS RATE_A]);
brchRateMva(:,3) = brchRateMva(:,3)/mpc_0.baseMVA;

%% Contingency sweep
% One contingency per generating unit: trip it, re-solve all three models,
% and count the violations each one reports.

numContingencies = size(mpc_0.gen,1);

% Per-contingency record
freq_SALT       = NaN(numContingencies,1);   % Hz
freq_enhPF      = NaN(numContingencies,1);   % Hz
outagedGenBus   = NaN(numContingencies,1);
Vload_SALT_all  = [];  Vload_enhPF_all = [];
Pgen_SALT_all   = [];  Pgen_enhPF_all  = [];

% Violation totals across the whole sweep
numVlimViolat_SALT  = 0;  numVlimViolat_enhPF  = 0;  numVlimViolat_PF  = 0;
numQlimViolat_SALT  = 0;  numQlimViolat_enhPF  = 0;  numQlimViolat_PF  = 0;
numLineViolat_SALT  = 0;  numLineViolat_enhPF  = 0;  numLineViolat_PF  = 0;
numNoConvergence_PF = 0;

for idxOutage = 1:numContingencies

    busNum_genLoss = mpc_0.gen(idxOutage,GEN_BUS);
    outagedGenBus(idxOutage) = busNum_genLoss;

    % ---------------------- Power flow side --------------------------
    % setGenContingency trips every unit at the bus, demotes it to PQ and,
    % if the slack was the unit lost, hands back the replacement slack so
    % SALT can be pointed at the same one.
    [mpc_contingency,newSlackBus] = setGenContingency(mpc_0,busNum_genLoss);

    [pfResult,~] = PF.solve(mpc_contingency,settings_nominal);
    bus_Vm_contingency_pf = pfResult.bus(:,VM);
    bus_Va_contingency_pf = deg2rad(pfResult.bus(:,VA));
    [I_gen_contingency_pf,I_load_contingency_pf,I_series_contingency_pf] = ...
        calcPfCurrents(pfResult);

    [enhPfResult,w_enhPF] = PF.solve(mpc_contingency,settings);
    bus_Vm_contingency_enhPf = enhPfResult.bus(:,VM);
    bus_Va_contingency_enhPf = deg2rad(enhPfResult.bus(:,VA));
    [I_gen_contingency_enhPf,I_load_contingency_enhPf, ...
        I_series_contingency_enhPf] = calcPfCurrents(enhPfResult);

    % The power flow keeps the tripped machine as a row carrying zero
    % current. SALT removes the device outright, so drop the row here too
    % and the two current lists stay comparable.
    for kList = 1:2
        if kList==1, Ig = I_gen_contingency_pf; else, Ig = I_gen_contingency_enhPf; end
        [~,idxDrop] = ismember(busNum_genLoss,Ig(:,1));
        assert(idxDrop>0 && Ig(idxDrop,2)==0)
        Ig(idxDrop,:) = [];
        if kList==1, I_gen_contingency_pf = Ig; else, I_gen_contingency_enhPf = Ig; end
    end

    % ------------------------- SALT side -----------------------------
    % SALT takes the outage as a BUS NUMBER and rebuilds from the base
    % case, so the device is genuinely absent rather than zeroed.
    salt = SALT(rawFile,deviceFile, ...
        maxDeltaT   = maxDeltaT, ...
        loadScale   = loadScale, ...
        lineOutage  = false, ...
        genOutage   = busNum_genLoss, ...
        loadOutage  = false, ...
        newSlackBus = newSlackBus, ...
        rerunPfInit = true);
    salt.caseName = char(caseName);

    scriptSimulate_SALT

    bus_VreIm_contingency_SALT = ...
        salt.x_v(salt.getRowIdxAbc_fromBusNums(salt.bus));
    bus_V_contingency_SALT = ...
        bus_VreIm_contingency_SALT(1:2:(2*salt.numBus)) + ...
     1j*bus_VreIm_contingency_SALT(1+(1:2:(2*salt.numBus)));
    bus_Vm_contingency_SALT = abs(bus_V_contingency_SALT);
    bus_Va_contingency_SALT = angle(bus_V_contingency_SALT);
    [I_gen_contingency_SALT,I_load_contingency_SALT, ...
        I_series_contingency_SALT] = calcSaltCurrents(salt);

    freq_SALT(idxOutage)  = salt.x_v(salt.getRowIdx_w)/(2*pi);
    freq_enhPF(idxOutage) = w_enhPF/(2*pi);

    % ------------------------ Derived quantities ---------------------
    % One call per model. The three models must agree on WHICH buses and
    % branches exist, or the comparison is not element for element.
    assert(isequal(I_gen_contingency_SALT(:,1),I_gen_contingency_enhPf(:,1)))
    assert(isequal(I_load_contingency_SALT(:,1),I_load_contingency_enhPf(:,1)))
    assert(isequal(I_series_contingency_SALT(:,1:2), ...
                   I_series_contingency_enhPf(:,1:2)))

    [Vload_salt,~,~,~,Pgen_salt] = postprocessData(bus, ...
        bus_Vm_contingency_SALT,bus_Va_contingency_SALT, ...
        I_gen_contingency_SALT,I_load_contingency_SALT, ...
        I_series_contingency_SALT);
    [Vload_enhPf,~,~,~,Pgen_enhPf] = postprocessData(bus, ...
        bus_Vm_contingency_enhPf,bus_Va_contingency_enhPf, ...
        I_gen_contingency_enhPf,I_load_contingency_enhPf, ...
        I_series_contingency_enhPf);

    Vload_SALT_all  = [Vload_SALT_all,  Vload_salt];  %#ok<AGROW>
    Vload_enhPF_all = [Vload_enhPF_all, Vload_enhPf]; %#ok<AGROW>
    Pgen_SALT_all   = [Pgen_SALT_all,   Pgen_salt];   %#ok<AGROW>
    Pgen_enhPF_all  = [Pgen_enhPF_all,  Pgen_enhPf];  %#ok<AGROW>

    % ------------------------- Violations ----------------------------
    % The outaged unit's limits go with it, so the limit lists are trimmed
    % the same way the current lists were.
    genQmin_contingency = genQmin;  genQmin_contingency(idxOutage,:) = [];
    genQmax_contingency = genQmax;  genQmax_contingency(idxOutage,:) = [];
    genPmin_contingency = genPmin;  genPmin_contingency(idxOutage,:) = [];
    genPmax_contingency = genPmax;  genPmax_contingency(idxOutage,:) = [];
    brchRateMva_contingency = brchRateMva;

    assert(converge_SS)
    [vlimViolat,busV_salt,~,~,qlimViolat,busQ_salt, ...
        ratedLineViolat,lineIdx_salt] = checkContingencyViolations( ...
        bus,bus_Vm_contingency_SALT,bus_Va_contingency_SALT, ...
        I_gen_contingency_SALT,I_load_contingency_SALT, ...
        I_series_contingency_SALT, ...
        genQmin_contingency,genQmax_contingency, ...
        genPmin_contingency,genPmax_contingency,busVmin,busVmax, ...
        brchRateMva_contingency);
    numVlimViolat_SALT = numVlimViolat_SALT + size(vlimViolat,1);
    numQlimViolat_SALT = numQlimViolat_SALT + size(qlimViolat,1);
    numLineViolat_SALT = numLineViolat_SALT + size(ratedLineViolat,1);

    assert(enhPfResult.success)
    [vlimViolat,busV_enhPf,~,~,qlimViolat,busQ_enhPf, ...
        ratedLineViolat,lineIdx_enhPf] = checkContingencyViolations( ...
        bus,bus_Vm_contingency_enhPf,bus_Va_contingency_enhPf, ...
        I_gen_contingency_enhPf,I_load_contingency_enhPf, ...
        I_series_contingency_enhPf, ...
        genQmin_contingency,genQmax_contingency, ...
        genPmin_contingency,genPmax_contingency,busVmin,busVmax, ...
        brchRateMva_contingency);
    numVlimViolat_enhPF = numVlimViolat_enhPF + size(vlimViolat,1);
    numQlimViolat_enhPF = numQlimViolat_enhPF + size(qlimViolat,1);
    numLineViolat_enhPF = numLineViolat_enhPF + size(ratedLineViolat,1);

    % The plain power flow is the one that can fail to converge, and a
    % contingency it cannot solve is counted separately rather than
    % silently scoring zero violations.
    if pfResult.success
        [vlimViolat,busV_pf,~,~,qlimViolat,busQ_pf, ...
            ratedLineViolat,lineIdx_pf] = checkContingencyViolations( ...
            bus,bus_Vm_contingency_pf,bus_Va_contingency_pf, ...
            I_gen_contingency_pf,I_load_contingency_pf, ...
            I_series_contingency_pf, ...
            genQmin_contingency,genQmax_contingency, ...
            genPmin_contingency,genPmax_contingency,busVmin,busVmax, ...
            brchRateMva_contingency);
        numVlimViolat_PF = numVlimViolat_PF + size(vlimViolat,1);
        numQlimViolat_PF = numQlimViolat_PF + size(qlimViolat,1);
        numLineViolat_PF = numLineViolat_PF + size(ratedLineViolat,1);
    else
        numNoConvergence_PF = numNoConvergence_PF + 1;
        busV_pf = []; busQ_pf = []; lineIdx_pf = zeros(0,2);
    end

    % One violation map per contingency
    plotContingencyViolations

end

%% Totals

fprintf('\n');
fprintf('================================================================\n');
fprintf(' Contingency analysis: %s\n',caseName);
fprintf(' %d single-generator outages\n',numContingencies);
fprintf('================================================================\n');
fprintf(' %-24s %8s %8s %8s\n','violations (total)','SALT','Enh-PF','PF');
fprintf(' %-24s %8d %8d %8d\n','bus voltage', ...
    numVlimViolat_SALT,numVlimViolat_enhPF,numVlimViolat_PF);
fprintf(' %-24s %8d %8d %8d\n','generator Q limit', ...
    numQlimViolat_SALT,numQlimViolat_enhPF,numQlimViolat_PF);
fprintf(' %-24s %8d %8d %8d\n','branch rating', ...
    numLineViolat_SALT,numLineViolat_enhPF,numLineViolat_PF);
fprintf(' %-24s %8d %8d %8d\n','ALL', ...
    numVlimViolat_SALT+numQlimViolat_SALT+numLineViolat_SALT, ...
    numVlimViolat_enhPF+numQlimViolat_enhPF+numLineViolat_enhPF, ...
    numVlimViolat_PF+numQlimViolat_PF+numLineViolat_PF);
fprintf('----------------------------------------------------------------\n');
fprintf(' %-24s %8s %8s %8d\n','contingencies unsolved','-','-', ...
    numNoConvergence_PF);
fprintf(' %-24s %8.4f %8.4f %8s\n','min frequency (Hz)', ...
    min(freq_SALT),min(freq_enhPF),'60.0000');
fprintf('================================================================\n');
