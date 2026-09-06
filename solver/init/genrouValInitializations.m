function [obj,genrou_VAbsPhA,genrou_IComplexPhA_busPU,THETA,genrou_efd,genrou_TorqueM] = ...
    genrouValInitializations(obj,bus_VComplexPhA,bus_Va)
% Set i (1:7), e (8:12), wr_delta, and theta_delta

[~,busIdx_genrou] = ismembertol(obj.genrou_bus,obj.bus);

% Calc voltage and current (line-to-line (LL), RMS, bus PU):
genrou_VComplexPhA_busPU = bus_VComplexPhA(busIdx_genrou);
genrou_SComplex_busPU = (obj.genrou_MW + 1j*obj.genrou_Mvar)/ ...
    obj.baseMVA;
genrou_IComplexPhA_busPU = ((genrou_SComplex_busPU')./ ...
    (genrou_VComplexPhA_busPU')).';
% Sanity check on current and power
assert(all(abs( ...
    abs(genrou_IComplexPhA_busPU) ...
    - ...
    sqrt(obj.genrou_MW.^2 + obj.genrou_Mvar.^2)/obj.baseMVA ...
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
% delta - Angle btwn dq-axis (rotor) volt and stator volt (bus angle)
% theta - Angle btwn rotor (dq-axis) and stator
genrou_SAbs = sqrt(obj.genrou_MW.^2 + obj.genrou_Mvar.^2);
PHI = sign(obj.genrou_MW./obj.genrou_Mvar - eps).*acos(obj.genrou_MW./genrou_SAbs);
Xq = -squeeze(obj.genrou_LMatrix(2,2,:));
Ra =  squeeze(obj.genrou_RMatrix(1,1,:));
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
Xd = -squeeze(obj.genrou_LMatrix(1,1,:));
Lad = squeeze(obj.genrou_LMatrix(1,4,:));
ifd = (eq + Xd.*id + Ra.*iq)./Lad;

% Controls - set as parameters in case does not include controls
% efd
Rfd = -squeeze(obj.genrou_RMatrix(4,4,:));
genrou_efd = Rfd.*ifd;
% TorqueM
phiq = squeeze(obj.genrou_LspMatrix(1,2,:)).*iq;
phid = squeeze(obj.genrou_LspMatrix(2,1,:)).*id ...
     + squeeze(obj.genrou_LspMatrix(2,4,:)).*ifd;
genrou_TorqueM = phid.*iq + phiq.*id;


% Final initialize
rowIdxGenrou = obj.getRowIdx_fromGenrouIndices(1:obj.numGenrou);
obj.x_v(rowIdxGenrou) = id; % id
obj.x_v(rowIdxGenrou+1) = iq; % iq
%obj.x_v(rowIdxGenrou+2) = 0; % i0
obj.x_v(rowIdxGenrou+3) = ifd; % ifd
%obj.x_v(rowIdxGenrou+4) = 0; % i1d
%obj.x_v(rowIdxGenrou+5) = 0; % i1q
%obj.x_v(rowIdxGenrou+6) = 0; % i2q
obj.x_v(rowIdxGenrou+7) = ed; % ed
obj.x_v(rowIdxGenrou+8) = eq; % eq
%obj.x_v(rowIdxGenrou+9) = 0; % e0
obj.x_v(rowIdxGenrou+10) = genrou_efd; % efd
obj.x_v(rowIdxGenrou+11) = obj.w0; % wr (assume nominal)
obj.x_v(rowIdxGenrou+12) = THETA; % theta
obj.x_v(rowIdxGenrou+16) = genrou_TorqueM;

% Set abc for current
for idxGenrou = 1:obj.numGenrou
    % Calculate iabc
    [obj.x_v(rowIdxGenrou(idxGenrou)+(13:15)),iabc_genrou,invT] = calcIabcFromIdq0( ...
        THETA(idxGenrou),id(idxGenrou),iq(idxGenrou), ...
        obj.genrouToBus_I(busIdx_genrou(idxGenrou),idxGenrou));

    % Sanity check on Iabc
    iabc_pf_genrou = real(genrou_IComplexPhA(idxGenrou)*[ ...
        1;(-0.5-1j*0.5*sqrt(3));(-0.5+1j*0.5*sqrt(3))]);
    assert(max(abs((iabc_genrou-iabc_pf_genrou)./iabc_pf_genrou))<1e-3)

    % Sanity check on Vabc
    vabc_pf_genrou = real(genrou_VComplexPhA(idxGenrou)*[ ...
        1;(-0.5-1j*0.5*sqrt(3));(-0.5+1j*0.5*sqrt(3))]);
    vabc_genrou = invT*[ed(idxGenrou);eq(idxGenrou);0];
    assert(max(abs((vabc_genrou-vabc_pf_genrou)./vabc_pf_genrou))<1e-3)

end

% Sanity check for relation between MW and MVar
[genrou_MW_calc,genrou_Mvar_calc] = calcPowerFromVoltCurrent( ...
    genrou_VComplexPhA_busPU,genrou_IComplexPhA_busPU,obj.baseMVA);
assert(all(abs(genrou_MW_calc-obj.genrou_MW)<obj.absTolCheck))
assert(all(abs(genrou_Mvar_calc-obj.genrou_Mvar)<obj.absTolCheck))
% As a P.U. note, this is equivalent

end
