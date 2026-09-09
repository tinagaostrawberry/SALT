function [yVals,iVals,Pval,Qval] = stampOfPfLoadStdMatrices_SS( ...
    ... % Current guess (rows are devices)
    Vr,Vi,w, ...
    ... % Device constants
    P,Q,IMpowerPercent,IMpowerPercentReact,dPdw,dQdw, ...
    ... % System parameters
    w0,useFreqDeviat)
%STAMPOFPFLOADSTDMATRICES_SS  Vectorised constant-power load stamp.
%
%   The standard load holds P and Q fixed against voltage, so the only
%   nonlinearity is the 1/|V|^2 in the injected current. The motor share of
%   the demand is still split out, because that is the only part the
%   frequency factor is allowed to act on.
%
%   Pval/Qval are the demands actually drawn at the current guess, handed
%   back for the result post-processing.

Vmag2 = Vr.^2 + Vi.^2;

P_ZIP = P.*(1 - IMpowerPercent);
P_IM  = P.*IMpowerPercent;
Q_ZIP = Q.*(1 - IMpowerPercentReact);
Q_IM  = Q.*IMpowerPercentReact;

if useFreqDeviat
    deltaW   = w - w0;
    freqFacP = 1 + dPdw.*deltaW;
    freqFacQ = 1 + dQdw.*deltaW;
else
    freqFacP = 1;
    freqFacQ = 1;
end
Pval = P_ZIP + P_IM.*freqFacP;
Qval = Q_ZIP + Q_IM.*freqFacQ;

% --- Real current equation ---------------------------------------------
IrNum      = Pval.*Vr + Qval.*Vi;
FIr        = IrNum ./ Vmag2;
dNumIr_dVr = Pval;
dFIr_dVr   = (Vmag2.*dNumIr_dVr - IrNum.*2.*Vr) ./ (Vmag2.^2);
dNumIr_dVi = Qval;
dFIr_dVi   = (Vmag2.*dNumIr_dVi - IrNum.*2.*Vi) ./ (Vmag2.^2);

% --- Imaginary current equation ----------------------------------------
IiNum      = Pval.*Vi - Qval.*Vr;
FIi        = IiNum ./ Vmag2;
dNumIi_dVr = -Qval;
dFIi_dVr   = (Vmag2.*dNumIi_dVr - IiNum.*2.*Vr) ./ (Vmag2.^2);
dNumIi_dVi = Pval;
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
