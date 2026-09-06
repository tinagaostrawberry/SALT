function [all_stampVal,I_v] = stampOfExcDc1Matrices_TD_flat( ...
    obj,all_stampVal,I_v,x_v,x_p)
% ---- EXCITER stamp, vectorized over all exciters (flat script form) -----
% Replaces the  for excIdx = 1:numExcDc1  loop in EMT.stampNonlinear_G_I and
% the two-layer call
%   obj.stampNonlinear_G_I_exciter -> stampOfExcMatrices_sparseIndices.
%
% MUST run after the generator block: it overwrites the generator's Efd row.
%
% LAYOUT (PSexc_(:,k) = absolute slot of local stamp k for every exciter):
%    1: 2  state 1 (Efd filter)   [-(Ke/2+Te/dT) , 1/2]
%    3: 4  state 2 (Vt sensor)    [1/2+Tr/dT , -1/2]
%    5: 7  state 3 (lead/lag)     1/2+[Tc Tb Tc]/dT
%    8: 9  state 4 (amplifier)    [-Ka/2 , 1/2+Ta/dT]
%   10:11  state 5 (rate fdbk)    [-(Kf1/dT) , 1/2+Tf1/dT]
%   12:18  "outside"  [dEfd/dx1 , dEfd/dEfd , dEfd/dwr , dVt/dvt , dVt/dvabc(1x3)]
%
% Everything here is scalar per device, so the whole block is plain
% elementwise arithmetic and is bit-identical to the loop.  The two
% 3-element dot products (|vabc|^2 and dVt/dvabc * vabc) are written as
% right-associated sums, which is what MATLAB's 3-term dot kernel produces.
% Function form (not a script): identical arithmetic, but MATLAB can
% JIT-compile it, which scripts running in a dynamic workspace cannot.

eN_ = obj.numExcDc1;

% ---- Index maps, whole population in one shot ---------------------------
eRW_  = obj.getRowIdx_fromExcDc1Indices((1:eN_).')-1;                        % N x 1
ePS_  = obj.getStampIdx_fromExcDc1Indices((1:eN_).')-1 + ...
    (1:obj.excDc1_numStampVals);
eROW_ = eRW_ + (1:obj.numExcDc1Eqns);
[~,eGi_] = ismembertol(obj.excDc1_bus,obj.genrou_bus);                          % gen of each exc
eGRW_ = obj.getRowIdx_fromGenrouIndices(eGi_(:))-1;
[~,~,eXEfd_,eXWr_] = getGenrouIndicesTD;
eVab_ = obj.getRowIdxAbc_fromBusNums(obj.excDc1_bus.');                      % 3 x N
[~,eBus_] = ismembertol(obj.excDc1_bus,obj.bus);

% ---- States -------------------------------------------------------------
e_x1   = x_v(eRW_+1);
e_efd  = x_v(eGRW_+eXEfd_);
e_vt   = x_v(eRW_+6);
e_vabc = x_v(eVab_);                                                      % 3 x N
e_wr   = x_v(eGRW_+eXWr_);
%
e_x1p  = x_p(eRW_+1);
e_x2p  = x_p(eRW_+2);
e_x3p  = x_p(eRW_+3);
e_x4p  = x_p(eRW_+4);
e_x5p  = x_p(eRW_+5);
e_vtp  = x_p(eRW_+6);

% ---- Constants ----------------------------------------------------------
eTr_=obj.excDc1_Tr(:);   eTa_=obj.excDc1_Ta(:);   eTc_=obj.excDc1_Tc(:);
eTb_=obj.excDc1_Tb(:);   eTe_=obj.excDc1_Te(:);   eTf1_=obj.excDc1_Tf1(:);
eKf1_=obj.excDc1_Kf1(:); eKa_=obj.excDc1_Ka(:);   eKe_=obj.excDc1_Ke(:);
eVref_=obj.excDc1_vref(:); eG2E_=obj.excDc1_genToExc(:);
eWs_ = obj.w0;   eDt_ = obj.deltaT;
eB2G_ = obj.busToGenrou_V(eBus_(:),eGi_(:)) + zeros(eN_,1);

e_I = zeros(eN_,obj.numExcDc1Eqns);

% ------------------------------ Linear -----------------------------------
% State 1:  -Ke*x1 - Te*dx1/dt + x4 = 0
all_stampVal(ePS_(:,1)) = -(eKe_/2+eTe_/eDt_);
all_stampVal(ePS_(:,2)) = 1/2;
e_I(:,1) = -(-(eKe_/2-eTe_/eDt_).*e_x1p + (1/2)*e_x4p);
% State 2:  x2 + Tr*dx2/dt - vt = 0
all_stampVal(ePS_(:,3)) = 1/2+eTr_/eDt_;
all_stampVal(ePS_(:,4)) = -1/2;
e_I(:,2) = -((1/2-eTr_/eDt_).*e_x2p-(1/2)*e_vtp);
% State 3:  lead/lag summing junction
all_stampVal(ePS_(:,5)) = 1/2+eTc_/eDt_;
all_stampVal(ePS_(:,6)) = 1/2+eTb_/eDt_;
all_stampVal(ePS_(:,7)) = 1/2+eTc_/eDt_;
e_I(:,3) = -((1/2-eTb_/eDt_).*e_x3p+(1/2-eTc_/eDt_).*(e_x2p+e_x5p)) + eVref_;
% State 4:  -Ka*x3 + x4 + Ta*dx4/dt = 0
all_stampVal(ePS_(:,8)) = -eKa_/2;
all_stampVal(ePS_(:,9)) = 1/2+eTa_/eDt_;
e_I(:,4) = -(-(eKa_/2).*e_x3p + (1/2-eTa_/eDt_).*e_x4p);
% State 5:  -Kf1*dx1/dt + x5 + Tf1*dx5/dt = 0
all_stampVal(ePS_(:,10)) = -(eKf1_/eDt_);
all_stampVal(ePS_(:,11)) = 1/2+eTf1_/eDt_;
e_I(:,5) = -((eKf1_/eDt_).*e_x1p + (1/2-eTf1_/eDt_).*e_x5p);

% ----------------------------- Nonlinear ---------------------------------
% Efd normalised by speed:  genToExc*genrou_efd - x1*wr_pu = 0
e_JEfd_x1  = -e_wr/eWs_;
e_JEfd_efd = eG2E_;
e_JEfd_wr  = -e_x1/eWs_;
e_FEfd     = e_efd.*eG2E_ - e_x1.*e_wr/eWs_;
e_IgenEfd  = e_JEfd_x1.*e_x1 + e_JEfd_efd.*e_efd + e_JEfd_wr.*e_wr + -e_FEfd;

% Vt:  -Vt + sqrt(2/3)*(va^2+vb^2+vc^2)^(1/2) = 0
e_vG = (eB2G_.'.*e_vabc).';                                               % N x 3, gen PU
e_sq = e_vG(:,1).*e_vG(:,1) + ...
    (e_vG(:,2).*e_vG(:,2) + e_vG(:,3).*e_vG(:,3));
e_JVt_dvt   = -1;
e_JVt_dvabc = sqrt(2/3)*( e_sq.^(-1/2) ).*e_vG.*eB2G_;                    % N x 3
e_FVt       = -e_vt + sqrt(2/3)*( e_sq.^(1/2) );
e_I(:,6) = e_JVt_dvt*e_vt + ...
    ( e_JVt_dvabc(:,1).*e_vG(:,1) + ...
     (e_JVt_dvabc(:,2).*e_vG(:,2) + e_JVt_dvabc(:,3).*e_vG(:,3)) ) + -e_FVt;

% ---- "Outside" index, concatenated at the end ---------------------------
all_stampVal(ePS_(:,12)) = e_JEfd_x1;
all_stampVal(ePS_(:,13)) = e_JEfd_efd;
all_stampVal(ePS_(:,14)) = e_JEfd_wr;
all_stampVal(ePS_(:,15)) = e_JVt_dvt;
all_stampVal(ePS_(:,16:18)) = e_JVt_dvabc;

% ---- Scatter the RHS ----------------------------------------------------
I_v(eROW_.')        = e_I.';
I_v(eGRW_+eXEfd_)   = e_IgenEfd;
end
