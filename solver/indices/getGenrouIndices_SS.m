function [idxIdq,idxIfd,idxEdq,idxEfd,idxWr,idxTheta,idxIphABC,idxTm] = ...
    getGenrouIndices_SS
%GETGENROUINDICES_SS  Accessor for the GENROU steady-state variable layout.
%   Reads the globals set by SALTINITINDICES, so there is a single source
%   of truth for the layout.
global idxGenrou_Idq idxGenrou_Ifd idxGenrou_Edq idxGenrou_Efd ...
       idxGenrou_Wr idxGenrou_Theta idxGenrou_IphABC idxGenrou_Tm %#ok<*GVMIS>
if isempty(idxGenrou_Idq), saltInitIndices; end
idxIdq = idxGenrou_Idq;      idxIfd    = idxGenrou_Ifd;
idxEdq = idxGenrou_Edq;      idxEfd    = idxGenrou_Efd;
idxWr  = idxGenrou_Wr;       idxTheta  = idxGenrou_Theta;
idxIphABC = idxGenrou_IphABC; idxTm    = idxGenrou_Tm;
end
