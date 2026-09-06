function [genrou_MW_noPU,genrou_Mvar_noPU] = calcPowerFromVoltCurrent(bus_v,bus_i,baseMVA)

% Should be complex
assert(all(~isreal([bus_v;bus_i])))

% Calculate angle
phi = angle(bus_v)-angle(bus_i);

% Calculate P and Q
% P = |V|*|I|*cos(phi)
% Q = |V|*|I|*sin(phi)
genrou_MW_noPU = abs(bus_v).*abs(bus_i).*cos(phi) * baseMVA;
genrou_Mvar_noPU = abs(bus_v).*abs(bus_i).*sin(phi) * baseMVA;

end