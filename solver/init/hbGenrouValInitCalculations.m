function [id,iq,ifd,ed,eq,efd,THETA,TorqueM,VAbsPhA,IComplexPhA_busPU] = ...
    hbGenrouValInitCalculations(obj,bus_VComplexPhA,bus_Va)
%HBGENROUVALINITCALCULATIONS  Steady-state GENROU states from the power flow.
%
%   Returns the machine's operating point as TIME-DOMAIN values; the caller
%   converts them to the unnormalized DFT convention. Everything sits at
%   DC, which is why there are no imaginary parts here.
%
%   Only the four windings that carry current at steady state are solved
%   for: id, iq and ifd. The dampers (i1d, i1q, i2q) and the zero sequence
%   (i0) are identically zero at a balanced single-frequency operating
%   point, because their flux is constant so no voltage is induced. HB
%   keeps them as unknowns anyway - once harmonics appear they are NOT
%   zero, and they are the dominant harmonic path through the machine - but
%   at initialisation they start from zero.
%
%   Angles follow Kundur section 3.6.2:
%     PHI    between stator voltage and current
%     DELTA  between the rotor q axis and the terminal voltage
%     THETA  between the rotor d axis and stator winding A

assert(all(abs(angle(bus_VComplexPhA) - bus_Va)<1e-9))
[~,busIdx] = ismembertol(obj.genrou_bus,obj.bus);
nG = obj.numGenrou;

% Terminal voltage and current, line-to-line RMS, bus per-unit
VComplexPhA_busPU = bus_VComplexPhA(busIdx);
SComplex_busPU = (obj.genrou_MW0 + 1j*obj.genrou_Mvar0)/obj.baseMVA;
IComplexPhA_busPU = ((SComplex_busPU')./(VComplexPhA_busPU')).';
assert(all(abs( abs(IComplexPhA_busPU) ...
    - sqrt(obj.genrou_MW0.^2 + obj.genrou_Mvar0.^2)/obj.baseMVA ...
    ./abs(VComplexPhA_busPU) )<1e-6))

% Into machine per-unit
IComplexPhA = IComplexPhA_busPU.*obj.busToGen_I(busIdx,(1:nG).');
IAbsPhA = abs(IComplexPhA);
VComplexPhA = VComplexPhA_busPU.*obj.busToGen_V(busIdx,(1:nG).');
VAbsPhA = abs(VComplexPhA);                 % Et, the terminal voltage

% Angles
SAbs = sqrt(obj.genrou_MW0.^2 + obj.genrou_Mvar0.^2);
PHI = sign(obj.genrou_Mvar0).*acos(obj.genrou_MW0./SAbs); % Not sure if sign is correct?
Xq = obj.genrou_Xq;
Ra = obj.genrou_Ra;
% NOTE: ParaEmt has sign of Q here, but then it does not meet
% HB constraints
DELTA = atan( ...
    (Xq.*IAbsPhA.*cos(PHI)-Ra.*IAbsPhA.*sin(PHI)) ...
    ./ ...
    (VAbsPhA + Ra.*IAbsPhA.*cos(PHI)+Xq.*IAbsPhA.*sin(PHI)) );
THETA = DELTA + bus_Va(busIdx) - pi/2;

% dq values
id = IAbsPhA.*sin(DELTA+PHI);
iq = IAbsPhA.*cos(DELTA+PHI);
ed = VAbsPhA.*sin(DELTA);
eq = VAbsPhA.*cos(DELTA);

% Field circuit
Xd  = obj.genrou_Xd;
Xad = obj.genrou_Lad;
ifd = (eq + Xd.*id + Ra.*iq)./Xad;
efd = obj.genrou_Rfd.*ifd;

% Mechanical torque. The speed-flux matrix is read straight out of the
% full [7 x 7 x nG] LspMatrix:
%   row 1 of Lsp*i is -phiq, row 2 is +phid
% so entry (1,2) is the iq coefficient and (2,1),(2,4) are the id and ifd
% coefficients. The row-1 entry already carries the sign the torque
% expression below expects, which is why "phiq" here is the negated flux.
% At initialisation i1d = i1q = i2q = 0, so their columns drop out.
Lsp = obj.genrou_LspMatrix;
phiq = reshape(Lsp(1,2,:),nG,1).*iq;
phid = reshape(Lsp(2,1,:),nG,1).*id + reshape(Lsp(2,4,:),nG,1).*ifd;
TorqueM = phid.*iq + phiq.*id;

% Cross-check: the terminal voltage and current must reproduce the
% dispatched power. This is the one place the whole per-unit chain can be
% caught going wrong.
[MW_calc,Mvar_calc] = calcPowerFromVoltCurrent( ...
    VComplexPhA_busPU,IComplexPhA_busPU,obj.baseMVA);
assert(all(abs(MW_calc  -obj.genrou_MW0  )<obj.absTolCheck))
assert(all(abs(Mvar_calc-obj.genrou_Mvar0)<obj.absTolCheck))

end
