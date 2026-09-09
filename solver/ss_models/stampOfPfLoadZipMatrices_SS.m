function [yVals,iVals,Pval,Qval] = stampOfPfLoadZipMatrices_SS( ...
    ... % Current guess (rows are devices)
    Vr,Vi,w, ...
    ... % ZIP shares of the static part
    Pz_ZIP,Pi_ZIP,Pp_ZIP, Qz_ZIP,Qi_ZIP,Qp_ZIP, ...
    ... % ZIP shares of the motor part (the only part that feels frequency)
    Pz_IM,Pi_IM,Pp_IM, Qz_IM,Qi_IM,Qp_IM, ...
    ... % Device constants
    dPdw,dQdw, ...
    ... % System parameters
    w0,useFreqDeviat)
%STAMPOFPFLOADZIPMATRICES_SS  Vectorised ZIP-load stamp for every load.
%
%   The load is split into a static ZIP polynomial and a motor ZIP
%   polynomial. Only the motor part is scaled by the frequency factor,
%   which is why the two sets of coefficients are carried separately all
%   the way into the Jacobian rather than being summed up front.
%
%   Unlike the constant-power load, P and Q here depend on |V|, so the
%   current derivatives pick up the extra dP/d|V| and dQ/d|V| terms.
%
%   Pval/Qval are the demands actually drawn at the current guess, handed
%   back for the result post-processing.

Vmag2 = Vr.^2 + Vi.^2;
Vmag  = sqrt(Vmag2);

P_ZIP = Pz_ZIP.*Vmag2 + Pi_ZIP.*Vmag + Pp_ZIP;
P_IM  = Pz_IM .*Vmag2 + Pi_IM .*Vmag + Pp_IM;
Q_ZIP = Qz_ZIP.*Vmag2 + Qi_ZIP.*Vmag + Qp_ZIP;
Q_IM  = Qz_IM .*Vmag2 + Qi_IM .*Vmag + Qp_IM;

dP_ZIP_dV = 2*Pz_ZIP.*Vmag + Pi_ZIP;
dP_IM_dV  = 2*Pz_IM .*Vmag + Pi_IM;
dQ_ZIP_dV = 2*Qz_ZIP.*Vmag + Qi_ZIP;
dQ_IM_dV  = 2*Qz_IM .*Vmag + Qi_IM;

if useFreqDeviat
    deltaW   = w - w0;
    freqFacP = 1 + dPdw.*deltaW;
    freqFacQ = 1 + dQdw.*deltaW;
else
    freqFacP = 1;
    freqFacQ = 1;
end
Pval  = P_ZIP + P_IM.*freqFacP;
Qval  = Q_ZIP + Q_IM.*freqFacQ;
dP_dV = dP_ZIP_dV + dP_IM_dV.*freqFacP;
dQ_dV = dQ_ZIP_dV + dQ_IM_dV.*freqFacQ;

% --- Real current equation ---------------------------------------------
IrNum      = Pval.*Vr + Qval.*Vi;
FIr        = IrNum ./ Vmag2;
dNumIr_dVr = Pval + Vr.*(dP_dV.*Vr./Vmag) + Vi.*(dQ_dV.*Vr./Vmag);
dFIr_dVr   = (Vmag2.*dNumIr_dVr - IrNum.*2.*Vr) ./ (Vmag2.^2);
dNumIr_dVi = Qval + Vr.*(dP_dV.*Vi./Vmag) + Vi.*(dQ_dV.*Vi./Vmag);
dFIr_dVi   = (Vmag2.*dNumIr_dVi - IrNum.*2.*Vi) ./ (Vmag2.^2);

% --- Imaginary current equation ----------------------------------------
IiNum      = Pval.*Vi - Qval.*Vr;
FIi        = IiNum ./ Vmag2;
dNumIi_dVr = -Qval + Vi.*(dP_dV.*Vr./Vmag) - Vr.*(dQ_dV.*Vr./Vmag);
dFIi_dVr   = (Vmag2.*dNumIi_dVr - IiNum.*2.*Vr) ./ (Vmag2.^2);
dNumIi_dVi = Pval + Vi.*(dP_dV.*Vi./Vmag) - Vr.*(dQ_dV.*Vi./Vmag);
dFIi_dVi   = (Vmag2.*dNumIi_dVi - IiNum.*2.*Vi) ./ (Vmag2.^2);

% --- RHS (linearised about the current operating point) ------------------
Irhs_Vr = -(FIr - dFIr_dVr.*Vr - dFIr_dVi.*Vi);
Irhs_Vi = -(FIi - dFIi_dVr.*Vr - dFIi_dVi.*Vi);

if useFreqDeviat
    dP_dw_val = P_IM.*dPdw;
    dQ_dw_val = Q_IM.*dQdw;
    dNumIr_dw = dP_dw_val.*Vr + dQ_dw_val.*Vi;
    dFIr_dw   = dNumIr_dw ./ Vmag2;
    dNumIi_dw = dP_dw_val.*Vi - dQ_dw_val.*Vr;
    dFIi_dw   = dNumIi_dw ./ Vmag2;

    Irhs_Vr = Irhs_Vr + dFIr_dw.*w;
    Irhs_Vi = Irhs_Vi + dFIi_dw.*w;

    % Order MUST match the index construction in PF.prepareStamping
    yVals = [dFIr_dVr; dFIr_dVi; dFIi_dVr; dFIi_dVi; ...
        dFIr_dw; dFIi_dw];
else
    yVals = [dFIr_dVr; dFIr_dVi; dFIi_dVr; dFIi_dVi];
end
iVals = [Irhs_Vr; Irhs_Vi];
end
