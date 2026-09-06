function saltInitIndices()
%SALTINITINDICES  Populate the STEADY-STATE index vectors as globals.
%
%   Call once - SALT does this from its constructor - and thereafter any
%   function can reach the index vectors with a single GLOBAL declaration
%   instead of re-calling a get*Indices function on every stamp.
%
%   Naming: idx<Model>_<state>, the position of each state inside that
%   device's block of SALT's unknown vector.
%
%   The time-domain layout is a different thing and lives in
%   EMTINITINDICES, which the EMT class calls from its own constructor. The
%   two differ because SALT carries every phasor as a (real,imag) pair
%   while EMT carries one instantaneous value per phase, so a solve that
%   never builds an EMT twin has no use for the TD vectors.
%
%   See also EMTINITINDICES.

% ---------------------------------------------------------------- GENROU
global idxGenrou_Idq idxGenrou_Ifd idxGenrou_Edq idxGenrou_Efd ...
       idxGenrou_Wr idxGenrou_Theta idxGenrou_IphABC idxGenrou_Tm       %#ok<*GVMIS>
idxGenrou_Idq    = 1:2;
idxGenrou_Ifd    = 3;
idxGenrou_Edq    = 4:5;
idxGenrou_Efd    = 6;
idxGenrou_Wr     = 7;
idxGenrou_Theta  = 8;
idxGenrou_IphABC = 9:14;
idxGenrou_Tm     = 15;

% --------------------------------------------------------------- EXCDC1
% 5 controller states + 1 terminal-voltage calculation
global idxExcDc1_x1 idxExcDc1_x2 idxExcDc1_x3 idxExcDc1_x4 idxExcDc1_x5 ...
       idxExcDc1_Vt
idxExcDc1_x1 = 1;  idxExcDc1_x2 = 2;  idxExcDc1_x3 = 3;
idxExcDc1_x4 = 4;  idxExcDc1_x5 = 5;  idxExcDc1_Vt = 6;

% --------------------------------------------------------------- TGOV1D
% 2 states + 1 mechanical-power calculation
global idxTgov1d_x1 idxTgov1d_x2 idxTgov1d_Pmech
idxTgov1d_x1 = 1;  idxTgov1d_x2 = 2;  idxTgov1d_Pmech = 3;

% ------------------------------------------------------------- GFM IBR
global idxIbrGFM_IphABC idxIbrGFM_Vg idxIbrGFM_Vo idxIbrGFM_If ...
       idxIbrGFM_Io idxIbrGFM_Vi idxIbrGFM_Gamma idxIbrGFM_IfRef ...
       idxIbrGFM_X idxIbrGFM_Wr idxIbrGFM_VoRefDroop_q idxIbrGFM_PqAvg ...
       idxIbrGFM_Theta idxIbrGFM_Vo_VR_avg idxIbrGFM_Vo_VR idxIbrGFM_VoRef
% Network interface
idxIbrGFM_IphABC        = 1:6;
% Passive network (LCL filter + coupling)
idxIbrGFM_Vg            = 7:8;
idxIbrGFM_Vo            = 9:10;
idxIbrGFM_If            = 11:12;
idxIbrGFM_Io            = 13:14;
idxIbrGFM_Vi            = 15:16;
% Current controller
idxIbrGFM_Gamma         = 17:18;
% Voltage controller
idxIbrGFM_IfRef         = 19:20;
idxIbrGFM_X             = 21:22;
% Droop
idxIbrGFM_Wr            = 23;
idxIbrGFM_VoRefDroop_q  = 24;
idxIbrGFM_PqAvg         = 25:26;
idxIbrGFM_Theta         = 27;
% Virtual resistor
idxIbrGFM_Vo_VR_avg     = 28:29;
idxIbrGFM_Vo_VR         = 30:31;
idxIbrGFM_VoRef         = 32:33;

% ------------------------------------------- COMPOSITE LOAD (IM + ZIP)
global idxLoadComp_Vs idxLoadComp_IphABC idxLoadComp_WSlip ...
       idxLoadComp_ThetaSys idxLoadComp_Is idxLoadComp_Ir ...
       idxLoadComp_Wr idxLoadComp_PQ idxLoadComp_Vm idxLoadComp_Izip
% Induction motor
idxLoadComp_Vs       = 1:2;
idxLoadComp_IphABC   = 3:8;
idxLoadComp_WSlip    = 9;
idxLoadComp_ThetaSys = 10;
idxLoadComp_Is       = 11:12;
idxLoadComp_Ir       = 13:14;
idxLoadComp_Wr       = 15;
% ZIP
idxLoadComp_PQ       = 16:17;
idxLoadComp_Vm       = 18;
idxLoadComp_Izip     = 19:20;

end
