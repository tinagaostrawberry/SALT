
function [idx,jdx,idx_dFParkV_dvabc,jdx_dFParkV_dvabc, ...
    idx_dFslip_dwSys,jdx_dFslip_dwSys,idx_dFsys_dwSys,jdx_dFsys_dwSys, ...
    idx_dFs_dwSys,jdx_dFs_dwSys,idx_dFswing_dwSys,jdx_dFswing_dwSys] = ...
    calcLoadCompIndices


% Indicies
[idxVs,idxIabc,idxWSlip,idxThetaSys,idxIs,idxIr,idxWr, ...
    idxPQ,idxVmPreFilt,idxVm,idxIzip] = getLoadCompIndicesTD;

% Prepare index storage
% (Start with reasonable estimate to avoid growth)
idxSizeEset = 1024;
idx = zeros(idxSizeEset,1);
jdx = zeros(idxSizeEset,1);
idxCount = 0;

% Helper to append a block (rows x cols) into triplet arrays
    function append_block(row,col)
        % Block is matrix of size length(rowIdx) x length(colIdx)
        % Append flattened triplets
        [rowExpand,colExpand] = ndgrid(row,col);
        m = numel(rowExpand);
        if idxCount + m > numel(idx)
            % grow by 2x
            newCap = max(numel(idx)*2, idxCount + m);
            idx(end+1:newCap) = 0;
            jdx(end+1:newCap) = 0;
        end
        idx(idxCount+1:idxCount+m) = rowExpand(:);
        jdx(idxCount+1:idxCount+m) = colExpand(:);
        idxCount = idxCount + m;
    end

% ------------------------- Induction motor -------------------------------

% Park
% vdq_s = T(theta_sys)*vabc
% => dFParkV_dVs_dThetaSys
%    dFParkV_dVabc
append_block(1:2,[idxVs,idxThetaSys])
[idx_dFParkV_dvabc,jdx_dFParkV_dvabc] = ndgrid(1:2,[0 0 0]);
idx_dFParkV_dvabc = idx_dFParkV_dvabc(:);
jdx_dFParkV_dvabc = jdx_dFParkV_dvabc(:);
% (idq_s+idq_ZIP) = T(theta_sys)*iabc
% 0 = [1 1 1]*iabc
% => dFParkI_dIs_dIzip_dThetaSys_dIabc
%    dFParkI0_dIabc
append_block(3:4,[idxIs,idxIzip,idxThetaSys,idxIabc])
append_block(5,idxIabc)

% Slip frequency
% wSlip = wSys - wr
% => dFslip_dwSlip_dwr
%    dFslip_dwSys
append_block(6,[idxWSlip,idxWr])
[idx_dFslip_dwSys,jdx_dFslip_dwSys] = ndgrid(6,0);

% System angle
% d(thetasys)/dt = wSys
% => dFsys_dThetaSys
%    dFsys_dwSys
append_block(7,idxThetaSys)
[idx_dFsys_dwSys,jdx_dFsys_dwSys] = ndgrid(7,0);

% Stator/rotor equations
% |vs| = [R]*|idq_s| + (d/dt)*[L]/w0*|idq_s| + [wMat]/w0*[L]*|idq_s|
% |0 |       |idq_r|                 |idq_r|                 |idq_r|
% s.t.
% R = |Rs_I2x2      |
%     |      Rr_I2x2|
%
% L = |Lss_I2x2 Lm_I2x2 |
%     |Lm_I2x2  Lrr_I2x2|
% => dFsr_dvs_disr_dwSlip
%    dFs_dwSys
append_block(8:11,[idxVs,idxIs,idxIr,idxWSlip])
[idx_dFs_dwSys,jdx_dFs_dwSys] = ndgrid(8:9,0);

% Swing equation
% (2*H/w0)*d(wr)/dt = ((isr^T)*L_Te*isr) - T0*(wr/wsys)^m
% => dFswing_dwr_disr
%    dFswing_dwSys = T0*wr*wSys^(-2);
append_block(12,[idxWr,idxIs,idxIr])
[idx_dFswing_dwSys,jdx_dFswing_dwSys] = ndgrid(12,0);

% ---------------------------- ZIP model ----------------------------------

% ZIP
% PQ = PQz*Vm^2 + PQi*Vm + PQp
% => dFZip_dPQ_dVm
append_block(13:14,[idxPQ,idxVm])

% Vm (pre-filt)
% VmPreFilt = (vd_s^2+vq_s^2)^(1/2)
% => dFVmPreFilt_dvmPreFilt_dvs
append_block(15,[idxVmPreFilt,idxVs])

% Vm (filt)
% Vm' = Tfilt*(VmPreFilt-Vm)
% => dFVm_dvm_dvmPreFilt
append_block(16,[idxVm,idxVmPreFilt])

% Izip
% idZip = (P*vd + Q*vq)/Vm^2
% iqZip = (P*vq - Q*vd)/Vm^2
% => dFIZip_dPQ_dvs_dVm_diZip
append_block(17:18,[idxPQ,idxVs,idxVm,idxIzip])


idx = idx(1:idxCount);
jdx = jdx(1:idxCount);

end