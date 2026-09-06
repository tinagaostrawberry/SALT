function [idx,jdx,idx_dFParkV_dvabc, jdx_dFParkV_dvabc] = calcIbrGFMIndices


% Indicies
[idxIabc,idxVg,idxVo,idxIf,idxIo,idxVi,idxGamma, ...
    idxIfRef,idxX,idxWr,idxVoRefDroop_q,idxPqAvg,idxTheta, ...
    idxVo_VR_avg,idxVo_VR,idxVoRef] = getIbrGFMIndicesTD;

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

% Park transform
% CURRENT:
% |io_dq| = |   T   | * |io_abc|
% | 0   |   | 1 1 1 |   |      |
% =>
% JibrGF(1:2,[idxIo idxTheta idxIabc]) = ...
%   [-eye(2) JT*iabc*busToIbrGFI T*busToIbrGFI];
% JibrGF(3,idxIabc) = 1;
append_block(1:2,[idxIo idxTheta idxIabc]);
append_block(3,idxIabc);
% VOLTAGE
% =>
% JibrGF(4:5,[idxVg idxTheta]) = [-eye(2) JT*vabc*busToIbrGFV];
append_block(4:5,[idxVg idxTheta]);
[idx_dFParkV_dvabc,jdx_dFParkV_dvabc] = ndgrid(4:5,[0 0 0]);
idx_dFParkV_dvabc = idx_dFParkV_dvabc(:);
jdx_dFParkV_dvabc = jdx_dFParkV_dvabc(:);

% Passive coupling
% io_dq' = (1/Lc)*(-rc*io_dq + vo_dq - vg_dq) + j*w*io_dq
append_block(6:7,[idxIo idxWr idxVo idxVg]);

% Passive capacitor
% vo_dq' = (1/Cf)*(if_dq - io_dq) + j*w*vo_dq + Rcap*(if_dq' - io_dq')
append_block(8:9,[idxVo idxWr idxIf idxIo]);

% Passive filter
% if_dq' = (1/Lf)*(-rf*if_dq + vi_dq - vo_dq) + j*w*if_dq
append_block(10:11,[idxIf idxWr idxVi idxVo]);

% Current controller
% 1. vi_dq = kc_i*gamma_dq + kc_p*d(gamma_dq)/dt - j*wr*Lf*if_dq + Gc*vo_dq
append_block(12:13,[idxVi idxGamma idxIf idxWr idxVo]);
% 2. d(gamma_dq)/dt = if_dq_ref - if_dq
append_block(14:15,[idxGamma idxIf idxIfRef]);

% Voltage controller
% 1. if_dq_ref = kv_i*x_dq + kv_p*d(x_dq)/dt - j*wr*Cf*vo_dq + Gv*io_dq
append_block(16:17,[idxIfRef idxX idxVo idxWr idxIo]);
% 2. d(x_dq)/dt = vo_dq_ref - vo_dq
append_block(18:19,[idxX idxVo idxVoRef]);

% Droop
% 1. wr = w + Mp*(pset - pavg)
append_block(20,[idxWr idxPqAvg(1)]);
% 2. vo_q_ref_droop = vo_q_set + Mq*(qset-qavg)
append_block(21,[idxVoRefDroop_q idxPqAvg(2)]);
% 3. pavg' = wmeas*(p-pavg) s.t. p = vg_dq.' * io_dq
%    NOTE: No factor of (3/2) is needed here. It would come from the
%    Park transform for peak, line-to-neutral, non-PU quantities; with
%    RMS, line-to-line, PU quantities it cancels.
append_block(22, [idxPqAvg(1) idxVg idxIo]);
% 4. qavg' = wmeas*(q-qavg) s.t. q = vg_dq.' * [0 -1; 1 0] * io_dq
append_block(23,[idxPqAvg(2) idxVg idxIo]);
% 5. d(theta)/dt = wr
append_block(24,[idxTheta idxWr]);

% Virtual resistor
% d/dt(Vo_VR_avg) = w_VR*(Rv*Io - Vo_VR_avg)
append_block(25:26,[idxVo_VR_avg idxIo]);
% Vo_VR = Vo_VR_avg - Rv*Io
append_block(27:28,[idxVo_VR idxVo_VR_avg idxIo]);
% Vo_ref = vo_ref_droop + Vo_VR
% =>
% JibrGF(29:30,[idxVoRef idxVo_VR]) = [-eye(2) eye(2)];
% JibrGF(30,idxVoRefDroop_q) = 1;
append_block(29:30,[idxVoRef idxVo_VR]);
append_block(30,idxVoRefDroop_q);

idx = idx(1:idxCount);
jdx = jdx(1:idxCount);

end
