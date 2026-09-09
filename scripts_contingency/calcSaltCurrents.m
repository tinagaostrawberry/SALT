function [I_gen,I_load,I_series] = calcSaltCurrents(salt)
%CALCSALTCURRENTS  Phase-A currents from a SALT solve.
%
%   [I_gen,I_load,I_series] = calcSaltCurrents_contingencyAnalysis_SS(salt)
%
%   salt      a CONVERGED SALT object.
%
%   I_gen     [busNum, I] per generating unit (GENROU then GFM IBR)
%   I_load    [busNum, I] per load (constant impedance then composite)
%   I_series  [busFrom, busTo, Iseries, Ifrom, Ito] per branch
%             (transmission lines then transformers)
%
%   All currents are the PHASE-A phasor in system per-unit, and every list
%   comes back sorted by bus number so it can be compared element by
%   element with the power-flow currents from
%   CALCPFCURRENTS_CONTINGENCYANALYSIS_SS.
%
%   The passive-branch currents are recomputed here from the converged
%   state with the SAME expressions STAMPLINEAR_Y_I uses, at the converged
%   system frequency - SALT always makes the linear network
%   frequency-dependent, so there is no nominal-frequency variant to pick
%   between.

w = salt.x_v(salt.getRowIdx_w);

%% Generating units

% ---- GENROU ----
% Terminal current is a state of the machine block; index 7 of the layout
% is the phase-A/B/C injection, of which entries 1-2 are phase A re/im.
[~,~,~,~,~,~,idxGenrou_IphABC] = getGenrouIndices_SS;
if salt.numGenrou > 0
    r = salt.getRowIdx_fromGenrouIndices((1:salt.numGenrou).')-1;
    I_genrou = [salt.genrou_bus(:), ...
        salt.x_v(r+idxGenrou_IphABC(1)) + 1j*salt.x_v(r+idxGenrou_IphABC(2))];
else
    I_genrou = zeros(0,2);
end

% ---- Grid-forming IBR ----
idxIbrGFM_IphABC = getIbrGFMIndices_SS;
if salt.numIbrGFM > 0
    r = salt.getRowIdx_fromIbrGFMIndices((1:salt.numIbrGFM).')-1;
    I_ibrGFM = [salt.ibrGFM_bus(:), ...
        salt.x_v(r+idxIbrGFM_IphABC(1)) + 1j*salt.x_v(r+idxIbrGFM_IphABC(2))];
else
    I_ibrGFM = zeros(0,2);
end

I_gen = [I_genrou;I_ibrGFM];

%% Loads

% ---- Constant impedance to ground ----
% Bus shunts and zero-real-power load buses. The admittance is the one the
% Jacobian used, so this is the same current the solver balanced.
nLd = numel(salt.load_bus);
if nLd > 0
    idx = salt.getRowIdxAbc_fromBusNums(salt.load_bus);
    % getRowIdxAbc_fromBusNums is phase-major ([phA(all); phB(all);
    % phC(all)]), so the gather runs bus-fastest within each phase:
    % reshape to (nLd x 3) and transpose, not straight to (3 x nLd).
    Vabc = reshape(complex(salt.x_v(idx(1:2:end)), ...
        salt.x_v(idx(2:2:end))),nLd,3).';
    Ym = salt.constLoadAdmittance(w);
    Iabc = Ym.*Vabc;
    I_loadZ = [salt.load_bus(:), Iabc(1,:).'];
else
    I_loadZ = zeros(0,2);
end

% ---- Composite load (induction motor + ZIP) ----
[~,idxLoadComp_IphABC] = getLoadCompIndices_SS;
if salt.numLoadComp > 0
    r = salt.getRowIdx_fromLoadCompIndices((1:salt.numLoadComp).')-1;
    I_loadComp = [salt.loadComp_bus(:), ...
        salt.x_v(r+idxLoadComp_IphABC(1)) + ...
        1j*salt.x_v(r+idxLoadComp_IphABC(2))];
else
    I_loadComp = zeros(0,2);
end

I_load = [I_loadZ;I_loadComp];

%% Series branches

% ---- Transmission lines (three-phase PI) ----
nTx = numel(salt.tx_busFrom);
if nTx > 0
    [VF,VT] = salt.branchVoltages(salt.tx_busFrom,salt.tx_busTo);
    Yb  = (1./(salt.tx_R + 1j*w*salt.tx_L)).';
    Ysh = (1j*w*salt.tx_C).';
    Ibr = Yb.*(VF-VT);
    Ifrom = Ibr + Ysh.*VF;          % leaves the "from" node
    Ito   = -Ibr + Ysh.*VT;         % leaves the "to" node
    I_tx = [salt.tx_busFrom(:),salt.tx_busTo(:), ...
        Ibr(1,:).',Ifrom(1,:).',Ito(1,:).'];
else
    I_tx = zeros(0,5);
end

% ---- Transformers ----
% Two kinds: a nonzero-resistance unit is stamped as an admittance, a
% lossless one carries its branch current as an MNA unknown. Both give the
% same [series, from, to] triple because a transformer has no shunt, so the
% current that enters one end leaves the other.
nXf = numel(salt.xfmr_busFrom);
if nXf > 0
    Iser = zeros(nXf,1);

    iA = salt.map.xfmr.idxAdm;
    if ~isempty(iA)
        [VF,VT] = salt.branchVoltages(salt.xfmr_busFrom(iA),salt.xfmr_busTo(iA));
        Ym = (1./(salt.xfmr_R(iA) + 1j*w*salt.xfmr_L(iA))).';
        Iabc = Ym.*(VF-VT);
        Iser(iA) = Iabc(1,:).';
    end

    iM = salt.map.xfmr.idxMna;
    if ~isempty(iM)
        Iabc = salt.mnaCurrents(iM);
        Iser(iM) = Iabc(1,:).';
    end

    I_xfmr = [salt.xfmr_busFrom(:),salt.xfmr_busTo(:),Iser,Iser,-Iser];
else
    I_xfmr = zeros(0,5);
end

I_series = [I_tx;I_xfmr];

%% Sort
% Sorted so the caller can compare against the power-flow lists row by row.

I_gen    = sortrows(I_gen,1);
I_load   = sortrows(I_load,1);
I_series = sortrows(I_series,[1 2]);

end
