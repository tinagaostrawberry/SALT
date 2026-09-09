function [Vload,Igen,SabsLine,PQload,Pgen,Vgen,Iload] = ...
    postprocessData_contingencyAnalysis_SS( ...
    bus,bus_Vm,bus_Va,I_gen,I_load,I_series)
%POSTPROCESSDATA_CONTINGENCYANALYSIS_SS  Derived quantities for ONE solve.
%
%   [Vload,Igen,SabsLine,PQload,Pgen,Vgen,Iload] = ...
%       postprocessData_contingencyAnalysis_SS( ...
%           bus,bus_Vm,bus_Va,I_gen,I_load,I_series)
%
%   bus       full bus-number list, in the order bus_Vm/bus_Va use
%   bus_Vm    bus voltage magnitude, pu
%   bus_Va    bus voltage angle, RADIANS
%   I_gen     [busNum, I] per generating unit
%   I_load    [busNum, I] per load
%   I_series  [busFrom, busTo, Iseries, Ifrom, Ito] per branch
%
%   Outputs, all magnitudes in per unit unless stated:
%     Vload     load-bus voltage magnitude
%     Igen      generator current magnitude
%     SabsLine  apparent power at both branch ends, [from; to]
%     PQload    load power, [P; Q] stacked
%     Pgen      generator real power (signed)
%     Vgen      generator-bus voltage magnitude
%     Iload     load current magnitude
%
%   This takes ONE model's results. Call it once per model (SALT, enhanced
%   power flow, plain power flow) and compare the returns; the earlier
%   two-model form made the caller pass one of the two sets twice whenever
%   it wanted a third model, and returned half its outputs unused.

%% Generator power

genBus = I_gen(:,1);
[~,idxGen] = ismember(genBus,bus);
assert(all(idxGen>0),'A generator bus is missing from the bus list')

VcomplexGen = bus_Vm(idxGen) .* exp(1j*bus_Va(idxGen));
IcomplexGen = I_gen(:,2);
Sgen = VcomplexGen .* conj(IcomplexGen);
Pgen = real(Sgen);

%% Load power

loadBus = I_load(:,1);
[~,idxLoad] = ismember(loadBus,bus);
assert(all(idxLoad>0),'A load bus is missing from the bus list')

VcomplexLoad = bus_Vm(idxLoad) .* exp(1j*bus_Va(idxLoad));
IcomplexLoad = I_load(:,2);
Sload = VcomplexLoad .* conj(IcomplexLoad);

%% Branch apparent power

fromBus = I_series(:,1);
toBus   = I_series(:,2);
[~,idxFrom] = ismember(fromBus,bus);
[~,idxTo]   = ismember(toBus,bus);
assert(all(idxFrom>0) && all(idxTo>0), ...
    'A branch terminal is missing from the bus list')

VcomplexFrom = bus_Vm(idxFrom) .* exp(1j*bus_Va(idxFrom));
VcomplexTo   = bus_Vm(idxTo)   .* exp(1j*bus_Va(idxTo));
% Columns 4 and 5 are the TOTAL currents at each end (series plus the line
% charging on that side), which is what a rating is written against.
SabsFrom = abs(VcomplexFrom .* conj(I_series(:,4)));
SabsTo   = abs(VcomplexTo   .* conj(I_series(:,5)));

%% Outputs

Vload    = abs(VcomplexLoad);
Igen     = abs(IcomplexGen);
SabsLine = [SabsFrom; SabsTo];
PQload   = [real(Sload); imag(Sload)];
Vgen     = abs(VcomplexGen);
Iload    = abs(IcomplexLoad);

end
