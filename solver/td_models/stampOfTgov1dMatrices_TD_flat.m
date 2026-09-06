function [all_stampVal,I_v] = stampOfTgov1dMatrices_TD_flat( ...
    obj,all_stampVal,I_v,x_v,x_p)
% ---- GOVERNOR stamp, vectorized over all governors (flat script form) ---
% Replaces the  for govIdx = 1:numTgov1d  loop in EMT.stampNonlinear_G_I and
% the two-layer call
%   obj.stampNonlinear_G_I_governor -> stampOfGovMatrices_sparseIndices.
%
% MUST run after the generator block: it overwrites the generator's Tm row.
%
% LAYOUT (PSgov_(:,k) = absolute slot of local stamp k for every governor):
%    1: 2  state 1  [1/2+T3/dT , -1/2-T2/dT]
%    3     state 2  -R/2 - R*T1/dT
%    4: 5  state 3  [1 , -1]
%    6:10  "outside" [dx2/dwr , dx3/dwr , dTm/dpmech , dTm/dwr , dTm/dTm]
%
% Every quantity is scalar per device, so this block is bit-identical to the
% loop.  The only product is the 3-term dTm/d[pmech;wr;Tm] * [pmech;wr;Tm]
% dot, written right-associated to match MATLAB's 3-term dot kernel.
% Function form (not a script): identical arithmetic, but MATLAB can
% JIT-compile it, which scripts running in a dynamic workspace cannot.

vN_ = obj.numTgov1d;

% ---- Index maps, whole population in one shot ---------------------------
vRW_  = obj.getRowIdx_fromTgov1dIndices((1:vN_).')-1;
vPS_  = obj.getStampIdx_fromTgov1dIndices((1:vN_).')-1 + ...
    (1:obj.tgov1d_numStampVals);
vROW_ = vRW_ + (1:obj.numTgov1dEqns);
[~,vGi_] = ismembertol(obj.tgov1d_bus,obj.genrou_bus);
vGRW_ = obj.getRowIdx_fromGenrouIndices(vGi_(:))-1;
[~,~,~,vXWr_,~,~,vXTm_] = getGenrouIndicesTD;

% ---- States -------------------------------------------------------------
v_pmech = x_v(vRW_+3);
v_Tm    = x_v(vGRW_+vXTm_);
v_wr    = x_v(vGRW_+vXWr_);
%
v_x1p   = x_p(vRW_+1);
v_x2p   = x_p(vRW_+2);
v_wrp   = x_p(vGRW_+vXWr_);

% ---- Constants ----------------------------------------------------------
vR_ =obj.tgov1d_R(:);  vT1_=obj.tgov1d_T1(:);  vT2_=obj.tgov1d_T2(:);
vT3_=obj.tgov1d_T3(:); vDt2_=obj.tgov1d_Dt(:); vPref_=obj.tgov1d_Pref(:);
vWs_ = obj.w0;  vDt_ = obj.deltaT;

v_I = zeros(vN_,obj.numTgov1dEqns);

% ------------------------------ Linear -----------------------------------
% State 1:  x1 + T3*dx1/dt - x2 - T2*dx2/dt = 0
all_stampVal(vPS_(:,1)) = 1/2+vT3_/vDt_;
all_stampVal(vPS_(:,2)) = -1/2-vT2_/vDt_;
v_I(:,1) = -((1/2-vT3_/vDt_).*v_x1p + (-1/2+vT2_/vDt_).*v_x2p);
% State 2:  pref - (wr/ws-1) - R*x2 - R*T1*dx2/dt = 0
all_stampVal(vPS_(:,3)) = -vR_/2 - vR_.*vT1_/vDt_;
v_JX2_dWr = -1/vWs_/2;
v_I(:,2) = -(vPref_ + (-vR_/2+vR_.*vT1_/vDt_).*v_x2p + -(1/vWs_/2)*v_wrp+1);
% State 3:  x1 - (wr/ws-1)*Dt - pmech = 0
all_stampVal(vPS_(:,4)) =  1;
all_stampVal(vPS_(:,5)) = -1;
v_JX3_dWr = -vDt2_/vWs_;
v_I(:,3) = -vDt2_;

% ----------------------------- Nonlinear ---------------------------------
% Tm:  pmech - (wr/ws)*Tm = 0
v_JTm_dPmech = ones(vN_,1);
v_JTm_dWr    = -v_Tm/vWs_;
v_JTm_dTm    = -v_wr/vWs_;
v_FTm        = v_pmech - (v_wr/vWs_).*v_Tm;
v_IgenTm = v_JTm_dPmech.*v_pmech + ...
    (v_JTm_dWr.*v_wr + v_JTm_dTm.*v_Tm) + -v_FTm;

% ---- "Outside" index, concatenated at the end ---------------------------
all_stampVal(vPS_(:,6))  = v_JX2_dWr + zeros(vN_,1);
all_stampVal(vPS_(:,7))  = v_JX3_dWr;
all_stampVal(vPS_(:,8))  = v_JTm_dPmech;
all_stampVal(vPS_(:,9))  = v_JTm_dWr;
all_stampVal(vPS_(:,10)) = v_JTm_dTm;

% ---- Scatter the RHS ----------------------------------------------------
I_v(vROW_.')      = v_I.';
I_v(vGRW_+vXTm_)  = v_IgenTm;
end
