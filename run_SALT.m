%% Steady-state After Last Transients (SALT)
% SALT is a steady-state solver whose solution is unified with the
% steady-state EMT response for transmission analysis. It computes its
% solution using physics-based EMT device models for single-harmonic,
% balanced systems of transmission networks. It solves for the complex bus
% voltages, internal device states, and system frequency, which varies from
% the 60 Hz nominal due to power mismatches.

%% Load data
% Network data comes from the power-flow case file; every dynamic device
% (GENROU, EXDC1, TGOV1, grid-forming IBR, ZIP load, induction motor) is
% described by a sheet in the workbook.

root = addpathSALT;

% Case table. Each row fixes the whole configuration of one test case, so
% the same numbers are used by run_SALT, by the reference capture in
% validation/reference and by the checks in validation.
%
%   scale     load perturbation, which DIVIDES the load power: below 1 is a
%             heavier load and therefore a lower system frequency.
%   genOut    bus number of the machine or inverter taken out, or false for
%             none. Every large case loses its THIRD generating unit, which
%             is what script_largeSystemTiming.m does.
%   verify    run the time-domain EMT check and plot it against SALT. That
%             costs about an hour per two periods at 2000 buses, so only
%             the small cases do it.
%
% The two 39-bus rows are the same network with different inverter
% penetration: 39bus keeps machines at buses 30, 31 and 39 (7 inverters),
% 39bus_highPenetrationIBR keeps only bus 30 (9 inverters).
%
% The six large systems are the cases in the timing driver
% script_largeSystemTiming.m.
caseName = "case300";

cases = struct( ...
    'name',   {"4bus","39bus","39bus_highPenetrationIBR","case118", ...
               "case300","case1354pegase","case3375wp","case2kbus", ...
               "case13659pegase"}, ...
    'net',    {'4bus.raw','case39_modified.m','case39_modified.m', ...
               'case118_modified.m','case300_modified.m', ...
               'case1354pegase_modified.m','case3375wp_modified.m', ...
               'case_ACTIVSg2000_modified.m','case13659pegase_modified.m'}, ...
    'dev',    {'4bus_SALT.xlsx','39bus_SALT.xlsx', ...
               '39bus_highPenetrationIBR_SALT.xlsx','case118_SALT.xlsx', ...
               'case300_SALT.xlsx','case1354pegase_SALT.xlsx', ...
               'case3375wp_SALT.xlsx','case2kbus_SALT.xlsx', ...
               'case13659pegase_SALT.xlsx'}, ...
    'plot',   {[], [30 31 34 39], [30 31 34 39], [1 12 25 100], ...
               [1 100 200 300], [], [], [7098 5262 1090 4192], []}, ...
    'scale',  {0.75, 0.97, 0.98, 1, 1, 1, 1, 1, 1}, ...
    'genOut', {false, false, 34, 6, 20, 221, 10080, 1009, 18}, ...
    'verify', {true, true, true, false, false, false, false, false, false});

idxCase = find(strcmp(string({cases.name}),caseName),1);
assert(~isempty(idxCase), ...
    'caseName must be one of: %s',strjoin(string({cases.name}),', '))
rawFile          = fullfile(root,'data',cases(idxCase).net);
deviceFile       = fullfile(root,'data',cases(idxCase).dev);
listBusToPlot    = cases(idxCase).plot;       % empty plots every bus
verifyWithEmtSim = cases(idxCase).verify;

%% Solver settings

% ----- Solver -----

% Newton-Raphson
totIters_NR = 50;
relErrTol_NR = 1e-4;
absErrTol_NR = 1e-8;

% Heuristics for convergence
include_clampLimiting = true;
clampLimit = 1;

% ----- Transient validation -----

% Periodicity check for the verification EMT run. Initialized from  SALT
% solution, so should be in steady-state, which is defined as waveforms
% being periodic for two periods within a relative error margin.
totIters_periodic = 50;
relErrTol_periodic = 1e-3;
absErrMagnitudePercentage_periodic = 1e-3;  % abs tol as a percentage of max

% Timestep for transient simulation, used after computing the steady-state
% solution via SALT to simulate and plot time-domain solution
% NOTE: Actual timestep size is largest number less than "maxDeltaT" that
% is a factor of the period of the waveform
% i.e. deltaT<=maxDeltaT, deltaT*numT=[1 period] s.t. numT is integer
maxDeltaT = 1e-5;

% Time-domain plotting
factor_downsample = 5;

%% Perturbation

% Contingencies. Each is false for none, or:
%   lineOutage  indices into the transmission-line list
%   genOutage   bus numbers whose machine or inverter is removed
%   loadOutage  bus numbers whose load is removed
lineOutage = false;
genOutage  = cases(idxCase).genOut;
loadOutage = false;

% Load perturbation scaling? (no perturbation is =1)
loadScale  = cases(idxCase).scale;

% Re-solve the power flow before initializing. See the README: with
% rerunPfInit true the base case is a converged operating point, so the
% SALT solution reproduces it at exactly nominal frequency.
rerunPfInit = true;

%% Build the solver
salt = SALT(rawFile,deviceFile, ...
    maxDeltaT   = maxDeltaT, ...
    loadScale   = loadScale, ...
    lineOutage  = lineOutage, ...
    genOutage   = genOutage, ...
    loadOutage  = loadOutage, ...
    rerunPfInit = rerunPfInit);

% Name the outputs after the TEST CASE, not the network file: the two
% 39-bus cases share case39_modified.m and would otherwise overwrite each
% other's plots and results file.
salt.caseName = char(caseName);

%% Solve, then verify
scriptSimulate_SALT
scriptPostprocess_SALT
