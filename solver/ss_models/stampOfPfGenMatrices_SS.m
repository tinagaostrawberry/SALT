function [yVals,iVals,P_real] = stampOfPfGenMatrices_SS( ...
    ... % Current guess (rows are devices)
    Vr,Vi,Qv,w, ...
    ... % Device constants
    P0,Vset,dPdw, ...
    ... % System parameters
    w0,useFreqDeviat)
%STAMPOFPFGENMATRICES_SS  Vectorised PV-generator stamp for every machine.
%
%   A PV generator contributes three equations at its bus: the real and
%   imaginary parts of the injected current, written against the unknown Q,
%   and the voltage-set-point constraint |V|^2 = Vset^2 that closes the
%   system for that extra Q unknown.
%
%   yVals is the flat list of Jacobian entries, iVals the flat list of RHS
%   entries. Both are ordered to match the triplet map PF.prepareStamping
%   builds, so they can be dropped straight into the value array without
%   any per-device indexing.
%
%   P_real is the dispatched active power actually seen by the network - it
%   differs from P0 only when frequency is a variable - and is handed back
%   so the result post-processing does not have to recompute it.

Vmag2 = Vr.^2 + Vi.^2;

if useFreqDeviat
    % Droop: the machine picks up load as the system slows down.
    deltaW = w - w0;
    P_real = P0 + dPdw.*deltaW;
else
    P_real = P0;
end

% --- Real current equation ---------------------------------------------
IrNum    = -P_real.*Vr + Qv.*Vi;
FIr      = IrNum ./ Vmag2;
dFIr_dQ  = Vi ./ Vmag2;
dFIr_dVr = (-Vmag2.*P_real - IrNum.*2.*Vr) ./ (Vmag2.^2);
dFIr_dVi = ( Vmag2.*Qv     - IrNum.*2.*Vi) ./ (Vmag2.^2);

% --- Imaginary current equation ----------------------------------------
IiNum    = -P_real.*Vi - Qv.*Vr;
FIi      = IiNum ./ Vmag2;
dFIi_dQ  = -Vr ./ Vmag2;
dFIi_dVr = (-Vmag2.*Qv     - IiNum.*2.*Vr) ./ (Vmag2.^2);
dFIi_dVi = (-Vmag2.*P_real - IiNum.*2.*Vi) ./ (Vmag2.^2);

% --- Voltage-set-point equation ----------------------------------------
FVg      = Vmag2 - Vset.^2;
dFVg_dVr = 2.*Vr;
dFVg_dVi = 2.*Vi;

% --- RHS (linearised about the current operating point) ------------------
% The solve is written as J(x_k)*x_{k+1} = J(x_k)*x_k - F(x_k), so the RHS
% is the residual with every Jacobian-times-state term added back.
Irhs_Vr = -(FIr - dFIr_dQ.*Qv - dFIr_dVr.*Vr - dFIr_dVi.*Vi);
Irhs_Vi = -(FIi - dFIi_dQ.*Qv - dFIi_dVr.*Vr - dFIi_dVi.*Vi);
Irhs_Q  = -(FVg - dFVg_dVr.*Vr - dFVg_dVi.*Vi);

if useFreqDeviat
    dFIr_dw = (-Vr.*dPdw) ./ Vmag2;
    dFIi_dw = (-Vi.*dPdw) ./ Vmag2;
    Irhs_Vr = Irhs_Vr + dFIr_dw.*w;
    Irhs_Vi = Irhs_Vi + dFIi_dw.*w;

    % Order MUST match the index construction in PF.prepareStamping
    yVals = [dFIr_dQ; dFIr_dVr; dFIr_dVi; ...
        dFIi_dQ; dFIi_dVr; dFIi_dVi; ...
        dFVg_dVr; dFVg_dVi; dFIr_dw; dFIi_dw];
else
    yVals = [dFIr_dQ; dFIr_dVr; dFIr_dVi; ...
        dFIi_dQ; dFIi_dVr; dFIi_dVi; ...
        dFVg_dVr; dFVg_dVi];
end
iVals = [Irhs_Vr; Irhs_Vi; Irhs_Q];
end
