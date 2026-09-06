function [idxVs,idxIphABC,idxWSlip,idxThetaSys,idxIs,idxIr,idxWr, ...
    idxPQ,idxVm,idxIzip] = getLoadCompIndices_SS
%GETLOADCOMPINDICES_SS  Accessor for the composite-load steady-state layout.
global idxLoadComp_Vs idxLoadComp_IphABC idxLoadComp_WSlip ...
       idxLoadComp_ThetaSys idxLoadComp_Is idxLoadComp_Ir ...
       idxLoadComp_Wr idxLoadComp_PQ idxLoadComp_Vm idxLoadComp_Izip %#ok<*GVMIS>
if isempty(idxLoadComp_Vs), saltInitIndices; end
idxVs = idxLoadComp_Vs; idxIphABC = idxLoadComp_IphABC;
idxWSlip = idxLoadComp_WSlip; idxThetaSys = idxLoadComp_ThetaSys;
idxIs = idxLoadComp_Is; idxIr = idxLoadComp_Ir; idxWr = idxLoadComp_Wr;
idxPQ = idxLoadComp_PQ; idxVm = idxLoadComp_Vm; idxIzip = idxLoadComp_Izip;
end
