function [idxI,idxE,idxEfd,idxWr,idxTheta,idxIabc,idxTm] = getGenrouIndicesTD
%GETGENROUINDICESTD  Accessor for the GENROU time-domain variable layout.
global idxGenrouTD_I idxGenrouTD_E idxGenrouTD_Efd idxGenrouTD_Wr ...
       idxGenrouTD_Theta idxGenrouTD_Iabc idxGenrouTD_Tm %#ok<*GVMIS>
if isempty(idxGenrouTD_I), emtInitIndices; end
idxI = idxGenrouTD_I;   idxE = idxGenrouTD_E;   idxEfd = idxGenrouTD_Efd;
idxWr = idxGenrouTD_Wr; idxTheta = idxGenrouTD_Theta;
idxIabc = idxGenrouTD_Iabc; idxTm = idxGenrouTD_Tm;
end
