function [vlimViolat,bus_voltViol, ...
    plimViolat,bus_P, ...
    qlimViolat,bus_Q, ...
    ratedLineViolat,busFromTo_ratedViolat] = ...
    checkContingencyViolations( ...
    bus, ...
    bus_Vm,bus_Va,I_gen,I_load,I_series, ...
    ...
    genQmin,genQmax,genPmin,genPmax,busVmin,busVmax, ...
    brchRateMva)
arguments
    bus
    bus_Vm
    bus_Va
    I_gen
    I_load %#ok<INUSA>
    I_series
    genQmin
    genQmax
    genPmin
    genPmax
    busVmin
    busVmax
    brchRateMva
end

%% Voltage violations?
assert(isequal(bus,busVmin(:,1)))
assert(isequal(bus,busVmax(:,1)))

% Separate min and max logicals to map the correct limit
idxV_min = bus_Vm < busVmin(:,2);
idxV_max = bus_Vm > busVmax(:,2);
idxV = logical(idxV_min + idxV_max);

bus_voltViol = bus(idxV);

% Construct violated value and limit pairs
v_lims = zeros(size(bus_Vm));
v_lims(idxV_min) = busVmin(idxV_min, 2);
v_lims(idxV_max) = busVmax(idxV_max, 2);
vlimViolat = [bus_Vm(idxV), v_lims(idxV)];

%% Power violations
genBus = I_gen(:,1);
assert(isequal(genQmin(:,1),genBus))
assert(isequal(genQmax(:,1),genBus))
assert(isequal(genPmin(:,1),genBus))
assert(isequal(genPmax(:,1),genBus))
[~,idxGen] = ismember(genBus,bus);

%
VcomplexGen = bus_Vm(idxGen) .* exp(1j*bus_Va(idxGen));
IcomplexGen = I_gen(:,2);
Sgen = VcomplexGen .* conj(IcomplexGen);
Pgen = real(Sgen);
Qgen = imag(Sgen);

% Separate logicals for P
idxP_min = Pgen < genPmin(:,2);
idxP_max = Pgen > genPmax(:,2);
idxP = logical(idxP_min + idxP_max);

% Separate logicals for Q
idxQ_min = Qgen < genQmin(:,2);
idxQ_max = Qgen > genQmax(:,2);
idxQ = logical(idxQ_min + idxQ_max);

bus_P = genBus(idxP);
bus_Q = genBus(idxQ);

% Construct P violated value and limit pairs
p_lims = zeros(size(Pgen));
p_lims(idxP_min) = genPmin(idxP_min, 2);
p_lims(idxP_max) = genPmax(idxP_max, 2);
plimViolat = [Pgen(idxP), p_lims(idxP)];

% Construct Q violated value and limit pairs
q_lims = zeros(size(Qgen));
q_lims(idxQ_min) = genQmin(idxQ_min, 2);
q_lims(idxQ_max) = genQmax(idxQ_max, 2);
qlimViolat = [Qgen(idxQ), q_lims(idxQ)];

%% Line violations
% Check all share same branch buses
fromBus = I_series(:,1);
toBus = I_series(:,2);
assert(isequal(fromBus,brchRateMva(:,1)))
assert(isequal(fromBus,I_series(:,1)))
assert(isequal(toBus,brchRateMva(:,2)))
assert(isequal(toBus,I_series(:,2)))

% Get indices
[~,idxFrom] = ismember(fromBus,bus);
[~,idxTo] = ismember(toBus,bus);

% Calculate complex voltages and currents
VcomplexFrom = bus_Vm(idxFrom) .* exp(1j*bus_Va(idxFrom));
VcomplexTo = bus_Vm(idxTo) .* exp(1j*bus_Va(idxTo));
IcomplexFrom = I_series(:,4);
IcomplexTo = I_series(:,5);

% Calculate complex power
SabsFrom = abs(VcomplexFrom .* conj(IcomplexFrom));
SabsTo = abs(VcomplexTo .* conj(IcomplexTo));

% Use the maximum apparent power between From and To for the violation value
SabsMax = max(SabsFrom, SabsTo);

% Find violation
idxRated = logical(SabsMax > brchRateMva(:,3));

% Get buses of violation
busFromTo_ratedViolat = [idxFrom(idxRated) idxTo(idxRated)];

% Construct line violated value and limit pairs
ratedLineViolat = [SabsMax(idxRated), brchRateMva(idxRated, 3)];

end