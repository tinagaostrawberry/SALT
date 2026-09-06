function [idxVs,idxIabc,idxWSlip,idxThetaSys,idxIs,idxIr,idxWr, ...
    idxPQ,idxVmPreFilt,idxVm,idxIzip] = getLoadCompIndicesTD
%GETLOADCOMPINDICESTD  Accessor for the composite-load time-domain layout.
global idxLoadCompTD_Vs idxLoadCompTD_Iabc idxLoadCompTD_WSlip ...
       idxLoadCompTD_ThetaSys idxLoadCompTD_Is idxLoadCompTD_Ir ...
       idxLoadCompTD_Wr idxLoadCompTD_PQ idxLoadCompTD_VmPreFilt ...
       idxLoadCompTD_Vm idxLoadCompTD_Izip %#ok<*GVMIS>
if isempty(idxLoadCompTD_Vs), emtInitIndices; end
idxVs = idxLoadCompTD_Vs; idxIabc = idxLoadCompTD_Iabc;
idxWSlip = idxLoadCompTD_WSlip; idxThetaSys = idxLoadCompTD_ThetaSys;
idxIs = idxLoadCompTD_Is; idxIr = idxLoadCompTD_Ir; idxWr = idxLoadCompTD_Wr;
idxPQ = idxLoadCompTD_PQ; idxVmPreFilt = idxLoadCompTD_VmPreFilt;
idxVm = idxLoadCompTD_Vm; idxIzip = idxLoadCompTD_Izip;
end
