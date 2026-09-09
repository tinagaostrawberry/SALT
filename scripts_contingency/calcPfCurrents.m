function [I_gen_pf,I_load_pf,I_series_pf] = calcPfCurrents(pfResult)
%CALCPFCURRENTS_CONTINGENCYANALYSIS_SS  Currents from a solved power flow.
%
%   [I_gen_pf,I_load_pf,I_series_pf] = ...
%       calcPfCurrents_contingencyAnalysis_SS(pfResult)
%
%   pfResult  a solved MATPOWER-format case (from PF.solve).
%
%   Same three lists, in the same units and the same sort order, as
%   CALCSALTCURRENTS_CONTINGENCYANALYSIS_SS, so the two can be compared
%   element by element.

define_constants;

%% Setup
% A bus-number -> row-index map, because bus numbers need not be 1..N.
nb = size(pfResult.bus, 1);
bus_ids = pfResult.bus(:, BUS_I);
b2i = sparse(bus_ids, 1, 1:nb);

% Complex bus voltage, per unit
V = pfResult.bus(:, VM) .* exp(1j * pfResult.bus(:, VA) * pi / 180);

%% Generator currents (injected into the bus)
g_bus_idx = b2i(pfResult.gen(:, GEN_BUS));
S_gen = (pfResult.gen(:, PG) + 1j * pfResult.gen(:, QG)) / pfResult.baseMVA;
I_gen_pf = [ ...
    pfResult.bus(g_bus_idx) ...
    conj(S_gen ./ V(g_bus_idx)) ...
    ];

%% Load currents (drawn from the bus)
S_load = (pfResult.bus(:, PD) + 1j * pfResult.bus(:, QD)) / pfResult.baseMVA;
I_load_all = conj(S_load ./ V);
load_mask = logical((pfResult.bus(:,PD)~=0) .* (pfResult.bus(:,QD)~=0));
assert(isequal(load_mask,(abs(I_load_all)>1e-6)))
I_load_pf = [ ...
    pfResult.bus(load_mask) ...
    I_load_all(load_mask) ...
    ];

%% Branch currents
f_idx = b2i(pfResult.branch(:, F_BUS));
t_idx = b2i(pfResult.branch(:, T_BUS));

R     = pfResult.branch(:, BR_R);
X     = pfResult.branch(:, BR_X);
Bc    = pfResult.branch(:, BR_B);       % total line charging susceptance
tap   = pfResult.branch(:, TAP);
shift = pfResult.branch(:, SHIFT);

% In MATPOWER a tap of 0 means "not a transformer", i.e. ratio 1.
tap(tap == 0) = 1;
tau = tap .* exp(1j * shift * pi / 180);

Ys  = 1 ./ (R + 1j * X);
Ysh = 1j * Bc / 2;                      % split equally between the ends

V_f = V(f_idx);
V_t = V(t_idx);

% Series current, referenced to the line side of the tap
I_series_val = ( (V_f ./ tau) - V_t ) .* Ys;

% Total current leaving the "from" node: Yff*Vf + Yft*Vt
I_from_total = ( (Ys + Ysh) ./ (abs(tau).^2) ) .* V_f ...
             - ( Ys ./ conj(tau) ) .* V_t;

% Total current leaving the "to" node: Ytf*Vf + Ytt*Vt
I_to_total = - ( Ys ./ tau ) .* V_f ...
           + ( Ys + Ysh ) .* V_t;

I_series_pf = [ ...
    pfResult.bus(f_idx) ...
    pfResult.bus(t_idx) ...
    I_series_val ...
    I_from_total ...
    I_to_total ...
    ];

%% Sort

I_gen_pf    = sortrows(I_gen_pf,1);
I_load_pf   = sortrows(I_load_pf,1);
I_series_pf = sortrows(I_series_pf,[1 2]);

end
