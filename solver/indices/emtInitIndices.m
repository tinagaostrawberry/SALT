function emtInitIndices()
%EMTINITINDICES  Populate the TIME-DOMAIN index vectors as globals.
%
%   Call once - EMT does this from its constructor - and thereafter any
%   function can reach the index vectors with a single GLOBAL declaration
%   instead of re-calling a get*IndicesTD function on every stamp.
%
%   Naming: idx<Model>TD_<state>, the position of each state inside that
%   device's block of EMT's unknown vector.
%
%   The steady-state layout is a different thing and lives in
%   SALTINITINDICES, which the SALT class calls from its own constructor.
%   The two differ because SALT carries every phasor as a (real,imag) pair
%   while EMT carries one instantaneous value per phase.
%
%   See also SALTINITINDICES.

% ---------------------------------------------------------------- GENROU
%#ok<*GVMIS>
global idxGenrouTD_I idxGenrouTD_E idxGenrouTD_Efd idxGenrouTD_Wr ...
       idxGenrouTD_Theta idxGenrouTD_Iabc idxGenrouTD_Tm
idxGenrouTD_I     = 1:7;
idxGenrouTD_E     = 8:10;
idxGenrouTD_Efd   = 11;
idxGenrouTD_Wr    = 12;
idxGenrouTD_Theta = 13;
idxGenrouTD_Iabc  = 14:16;
idxGenrouTD_Tm    = 17;

% ------------------------------------------------------------- GFM IBR
global idxIbrGFMTD_Iabc idxIbrGFMTD_Vg idxIbrGFMTD_Vo idxIbrGFMTD_If ...
       idxIbrGFMTD_Io idxIbrGFMTD_Vi idxIbrGFMTD_Gamma idxIbrGFMTD_IfRef ...
       idxIbrGFMTD_X idxIbrGFMTD_Wr idxIbrGFMTD_VoRefDroop_q ...
       idxIbrGFMTD_PqAvg idxIbrGFMTD_Theta idxIbrGFMTD_Vo_VR_avg ...
       idxIbrGFMTD_Vo_VR idxIbrGFMTD_VoRef
idxIbrGFMTD_Iabc          = 1:3;
idxIbrGFMTD_Vg            = 4:5;
idxIbrGFMTD_Vo            = 6:7;
idxIbrGFMTD_If            = 8:9;
idxIbrGFMTD_Io            = 10:11;
idxIbrGFMTD_Vi            = 12:13;
idxIbrGFMTD_Gamma         = 14:15;
idxIbrGFMTD_IfRef         = 16:17;
idxIbrGFMTD_X             = 18:19;
idxIbrGFMTD_Wr            = 20;
idxIbrGFMTD_VoRefDroop_q  = 21;
idxIbrGFMTD_PqAvg         = 22:23;
idxIbrGFMTD_Theta         = 24;
idxIbrGFMTD_Vo_VR_avg     = 25:26;
idxIbrGFMTD_Vo_VR         = 27:28;
idxIbrGFMTD_VoRef         = 29:30;

% ------------------------------------------- COMPOSITE LOAD (IM + ZIP)
global idxLoadCompTD_Vs idxLoadCompTD_Iabc idxLoadCompTD_WSlip ...
       idxLoadCompTD_ThetaSys idxLoadCompTD_Is idxLoadCompTD_Ir ...
       idxLoadCompTD_Wr idxLoadCompTD_PQ idxLoadCompTD_VmPreFilt ...
       idxLoadCompTD_Vm idxLoadCompTD_Izip
idxLoadCompTD_Vs        = 1:2;
idxLoadCompTD_Iabc      = 3:5;
idxLoadCompTD_WSlip     = 6;
idxLoadCompTD_ThetaSys  = 7;
idxLoadCompTD_Is        = 8:9;
idxLoadCompTD_Ir        = 10:11;
idxLoadCompTD_Wr        = 12;
idxLoadCompTD_PQ        = 13:14;
idxLoadCompTD_VmPreFilt = 15;
idxLoadCompTD_Vm        = 16;
idxLoadCompTD_Izip      = 17:18;

end
