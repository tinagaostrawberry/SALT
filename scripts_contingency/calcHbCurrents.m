function [I_gen,I_load,I_series] = calcHbCurrents(hb)
%CALCHBCURRENTS  Fundamental phase-A currents from an HB solve.
%
%   [I_gen,I_load,I_series] = calcHbCurrents(hb)
%
%   hb        a CONVERGED HB object.
%
%   I_gen     [busNum, I] per generating unit (GENROU then GFM IBR)
%   I_load    [busNum, I] per load (constant impedance then composite)
%   I_series  [busFrom, busTo, Iseries, Ifrom, Ito] per branch
%             (transmission lines then transformers)
%
%   The HB analogue of CALCSALTCURRENTS, element for element, so the two can
%   be compared row by row and either can be handed to
%   CHECKCONTINGENCYVIOLATIONS. Two things differ:
%
%   1. WHICH HARMONIC. SALT has one phasor per quantity; HB has a
%      coefficient per harmonic. The violation limits - bus voltage bands,
%      generator Q, branch MVA ratings - are all written against the
%      FUNDAMENTAL, so k = 1 is what is extracted here. Harmonic content, if
%      the run has any, is deliberately not folded in: doing so would
%      change what the limit means without saying so.
%
%   2. SCALING. x_v holds UNNORMALIZED DFT coefficients under
%      x(n) = (1/N)*sum_k X(k)*exp(1j*2*pi*k*n/N), and a real signal splits
%      a non-DC component across a conjugate pair, so the phasor is 2/N
%      times the stored coefficient. Everything below is converted once,
%      through SCALE, and then the arithmetic is identical to SALT's.
%
%   The passive-branch currents are recomputed from the converged state with
%   the same expressions STAMPLINEAR_Y_I uses, at the converged fundamental.
%
%   See also CALCSALTCURRENTS, CALCPFCURRENTS.

global idxHbGenrou_IphABC idxHbIbrGFM_IphABC ...
       idxHbLoadComp_IphABC                                    %#ok<GVMIS>

K_FUND = 1;                 % the harmonic the limits are written against
w = hb.x_v(hb.getRowIdx_w0);
wk = K_FUND*w;
% Coefficient -> phasor. DC would need N rather than N/2, but nothing here
% is a DC quantity.
scale = 2/hb.N;

% One quantity's phase-A phasor, for a whole device fleet
    function c = devPhA(rowBaseFcn,nDev,locIphABC)
        r = rowBaseFcn((1:nDev).',K_FUND)-1;
        c = scale*(hb.x_v(r+locIphABC(1)) + 1j*hb.x_v(r+locIphABC(2)));
    end

%% Generating units

% ---- GENROU ----
if hb.numGenrou > 0
    I_genrou = [hb.genrou_bus(:), ...
        devPhA(@(i,k)hb.getRowIdx_fromGenrouIndices(i,k), ...
        hb.numGenrou,idxHbGenrou_IphABC)];
else
    I_genrou = zeros(0,2);
end

% ---- Grid-forming IBR ----
if hb.numIbrGFM > 0
    I_ibrGFM = [hb.ibrGFM_bus(:), ...
        devPhA(@(i,k)hb.getRowIdx_fromIbrGFMIndices(i,k), ...
        hb.numIbrGFM,idxHbIbrGFM_IphABC)];
else
    I_ibrGFM = zeros(0,2);
end

I_gen = [I_genrou;I_ibrGFM];

%% Loads

% ---- Constant impedance to ground ----
nLd = numel(hb.load_bus);
if nLd > 0
    idx = hb.getRowIdxAbc_fromBusNums(hb.load_bus,K_FUND);
    % getRowIdxAbc_fromBusNums is phase-major ([phA(all); phB(all);
    % phC(all)]), so the gather runs bus-fastest within each phase:
    % reshape to (nLd x 3) and transpose, not straight to (3 x nLd).
    Vabc = scale*reshape(complex(hb.x_v(idx(1:2:end)), ...
        hb.x_v(idx(2:2:end))),nLd,3).';
    Ym = hb.constLoadAdmittance(wk,K_FUND);
    Iabc = Ym.*Vabc;
    I_loadZ = [hb.load_bus(:), Iabc(1,:).'];
else
    I_loadZ = zeros(0,2);
end

% ---- Composite load (induction motor + ZIP) ----
if hb.numLoadComp > 0
    I_loadComp = [hb.loadComp_bus(:), ...
        devPhA(@(i,k)hb.getRowIdx_fromLoadCompIndices(i,k), ...
        hb.numLoadComp,idxHbLoadComp_IphABC)];
else
    I_loadComp = zeros(0,2);
end

I_load = [I_loadZ;I_loadComp];

%% Series branches

% ---- Transmission lines (three-phase PI) ----
nTx = numel(hb.tx_busFrom);
if nTx > 0
    [VF,VT] = hb.branchVoltages(hb.tx_busFrom,hb.tx_busTo,K_FUND);
    VF = scale*VF;  VT = scale*VT;
    Yb  = (1./(hb.tx_R + 1j*wk*hb.tx_L)).';
    Ysh = (1j*wk*hb.tx_C).';
    Ibr = Yb.*(VF-VT);
    Ifrom = Ibr + Ysh.*VF;          % leaves the "from" node
    Ito   = -Ibr + Ysh.*VT;         % leaves the "to" node
    I_tx = [hb.tx_busFrom(:),hb.tx_busTo(:), ...
        Ibr(1,:).',Ifrom(1,:).',Ito(1,:).'];
else
    I_tx = zeros(0,5);
end

% ---- Transformers ----
% Two kinds: a nonzero-resistance unit is stamped as an admittance, a
% lossless one carries its branch current as an MNA unknown. Both give the
% same [series, from, to] triple because a transformer has no shunt, so the
% current that enters one end leaves the other.
nXf = numel(hb.xfmr_busFrom);
if nXf > 0
    Iser = zeros(nXf,1);

    iA = hb.map.xfmr.idxAdm;
    if ~isempty(iA)
        [VF,VT] = hb.branchVoltages(hb.xfmr_busFrom(iA), ...
            hb.xfmr_busTo(iA),K_FUND);
        Ym = (1./(hb.xfmr_R(iA) + 1j*wk*hb.xfmr_L(iA))).';
        Iabc = Ym.*scale.*(VF-VT);
        Iser(iA) = Iabc(1,:).';
    end

    iM = hb.map.xfmr.idxMna;
    if ~isempty(iM)
        Iabc = scale*hb.mnaCurrents(iM,K_FUND);
        Iser(iM) = Iabc(1,:).';
    end

    I_xfmr = [hb.xfmr_busFrom(:),hb.xfmr_busTo(:),Iser,Iser,-Iser];
else
    I_xfmr = zeros(0,5);
end

I_series = [I_tx;I_xfmr];

%% Sort
% Sorted so the caller can compare against the SALT and power-flow lists
% row by row.

I_gen    = sortrows(I_gen,1);
I_load   = sortrows(I_load,1);
I_series = sortrows(I_series,[1 2]);

end
