function [id,iq,ifd,ed,eq,genrou_efd,THETA,genrou_TorqueM, ...
    genrou_VAbsPhA,genrou_VComplexPhA_busPU,genrou_IComplexPhA_busPU] = ...
    genrouValInitCalculations_assumeW0(obj,bus_VComplexPhA,bus_Va)

assert(all(abs(angle(bus_VComplexPhA) - bus_Va)<1e-9))
[~,busIdx_genrou] = ismembertol(obj.genrou_bus,obj.bus);

% Calc voltage and current (line-to-line (LL), RMS, bus PU):
genrou_VComplexPhA_busPU = bus_VComplexPhA(busIdx_genrou);
genrou_SComplex_busPU = (obj.genrou_MW0 + 1j*obj.genrou_Mvar0)/ ...
    obj.baseMVA;
genrou_IComplexPhA_busPU = ((genrou_SComplex_busPU')./ ...
    (genrou_VComplexPhA_busPU')).';
% Sanity check on current and power
assert(all(abs( ...
    abs(genrou_IComplexPhA_busPU) ...
    - ...
    sqrt(obj.genrou_MW0.^2 + obj.genrou_Mvar0.^2)/obj.baseMVA ...
    ./abs(genrou_VComplexPhA_busPU) ...
    )<1e-6))

% Convert current to from bus to gen PU
genrou_IComplexPhA = genrou_IComplexPhA_busPU.* ...
    obj.busToGenrou_I(busIdx_genrou,1:obj.numGenrou);
genrou_IAbsPhA = abs(genrou_IComplexPhA);
genrou_VComplexPhA = genrou_VComplexPhA_busPU.* ...
    obj.busToGenrou_V(busIdx_genrou,1:obj.numGenrou);
genrou_VAbsPhA = abs(genrou_VComplexPhA); % Et (terminal voltage)

% Calc angles:
% phi   - Angle btwn stator volt and current
% delta - Angle btwn rotor q-axis and stator volt Et
% theta - Angle btwn rotor d-axis and stator winding A (changes wrt time)
genrou_SAbs = sqrt(obj.genrou_MW0.^2 + obj.genrou_Mvar0.^2);
PHI = sign(obj.genrou_Mvar0).*acos(obj.genrou_MW0./genrou_SAbs);
Xq = obj.genrou_Xq;
Ra = obj.genrou_Ra;
DELTA = atan( ...
    (Xq.*genrou_IAbsPhA.*cos(PHI)-Ra.*genrou_IAbsPhA.*sin(PHI)) ...
    ./ ...
    (genrou_VAbsPhA + ...
    Ra.*genrou_IAbsPhA.*cos(PHI)+Xq.*genrou_IAbsPhA.*sin(PHI)) ...
    );
THETA = DELTA + bus_Va(busIdx_genrou) - pi/2;

% Get DQ0 values
id = genrou_IAbsPhA.*sin(DELTA+PHI);
iq = genrou_IAbsPhA.*cos(DELTA+PHI);
ed = genrou_VAbsPhA.*sin(DELTA);
eq = genrou_VAbsPhA.*cos(DELTA);

% Get field values
Xd = obj.genrou_Xd;
Xad = obj.genrou_Lad;
ifd = (eq + Xd.*id + Ra.*iq)./Xad;

% Controls - set as parameters in case does not include controls
% efd
Rfd = obj.genrou_Rfd;
genrou_efd = Rfd.*ifd;
% TorqueM
phiq = squeeze(obj.genrou_LspMatrix(1,2,:)).*iq;
phid = squeeze(obj.genrou_LspMatrix(2,1,:)).*id ...
     + squeeze(obj.genrou_LspMatrix(2,4,:)).*ifd;
genrou_TorqueM = phid.*iq + phiq.*id;
% Another way to calculate torque
end