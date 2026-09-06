function [idxIphABC,idxVg,idxVo,idxIf,idxIo,idxVi,idxGamma, ...
    idxIfRef,idxX,idxWr,idxVoRefDroop_q,idxPqAvg,idxTheta, ...
    idxVo_VR_avg,idxVo_VR,idxVoRef] = getIbrGFMIndices_SS
%GETIBRGFMINDICES_SS  Accessor for the GFM IBR steady-state layout.
global idxIbrGFM_IphABC idxIbrGFM_Vg idxIbrGFM_Vo idxIbrGFM_If ...
       idxIbrGFM_Io idxIbrGFM_Vi idxIbrGFM_Gamma idxIbrGFM_IfRef ...
       idxIbrGFM_X idxIbrGFM_Wr idxIbrGFM_VoRefDroop_q idxIbrGFM_PqAvg ...
       idxIbrGFM_Theta idxIbrGFM_Vo_VR_avg idxIbrGFM_Vo_VR idxIbrGFM_VoRef %#ok<*GVMIS>
if isempty(idxIbrGFM_IphABC), saltInitIndices; end
idxIphABC = idxIbrGFM_IphABC; idxVg = idxIbrGFM_Vg; idxVo = idxIbrGFM_Vo;
idxIf = idxIbrGFM_If; idxIo = idxIbrGFM_Io; idxVi = idxIbrGFM_Vi;
idxGamma = idxIbrGFM_Gamma; idxIfRef = idxIbrGFM_IfRef; idxX = idxIbrGFM_X;
idxWr = idxIbrGFM_Wr; idxVoRefDroop_q = idxIbrGFM_VoRefDroop_q;
idxPqAvg = idxIbrGFM_PqAvg; idxTheta = idxIbrGFM_Theta;
idxVo_VR_avg = idxIbrGFM_Vo_VR_avg; idxVo_VR = idxIbrGFM_Vo_VR;
idxVoRef = idxIbrGFM_VoRef;
end
