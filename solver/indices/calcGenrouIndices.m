function [idx,jdx,idx_dFParkV_dvabc,jdx_dFParkV_dvabc] = calcGenrouIndices


% Indicies
[idxI,idxE,idxEfd,idxWr,idxTheta,idxIabc,idxTm] = getGenrouIndicesTD;

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

% Stator
% e = -[R]*i + [L]*(d/dt)i + [Lsp]*i*wr
% =>
% Jgen(idxI,[idxI idxWr]) = [-R+2*L./deltaT_PU+Lsp*wr/ws, Lsp*i/ws];
append_block(idxI,[idxI idxWr]);
% Jgen(1:4,[idxE idxEfd]) = -eye(4);
append_block(1:4,[idxE idxEfd]);

% Park for voltage
% vabc = [T(theta)^-1]*vdq0
% s.t.: (d/dt)theta = wr = (d/dt)d + ws => theta = d + ws*t
% =>
% Jgen(idxE,[idxE idxTheta]) = [invT, JinvT*e(1:3)];
append_block(idxE,[idxE idxTheta]);
[idx_dFParkV_dvabc,jdx_dFParkV_dvabc] = ndgrid(idxE,[0 0 0]);
idx_dFParkV_dvabc = idx_dFParkV_dvabc(:);
jdx_dFParkV_dvabc = jdx_dFParkV_dvabc(:);

% Field voltage control
% Depends on exc control, so leave row 11 of Jacobian empty for now!

% Swing
% (d/dt)(deltaWr_pu) = (1/(2*H))*(Tm - Te - KD*deltaWr)
% s.t.: Te = -phiq*id + phid*iq = i^T*([-phiq phid | 0 ])
%          = i^T*Lsp*i
% deltaWr_pu = (1/ws)*(wr - ws) = (1/ws)*((d/t)theta - ws)
% =>
% Jgen(idxWr,[idxI idxWr idxTm]) = [(Lsp+Lsp.')*i (4*H/deltaT+KD)/ws -1];
append_block(idxWr,[idxI idxWr idxTm]);
% =>
% Jgen(idxTheta,[idxWr idxTheta]) = [-deltaT 2];
append_block(idxTheta,[idxWr idxTheta]);

% Park for current
% iabc = [T(theta)^-1]*idq0
% =>
% Jgen(idxIabc,[(1:3) idxTheta idxIabc]) = ...
%   [invT JinvT*i(1:3) -eye(3)*busToGenI];
append_block(idxIabc,[(1:3) idxTheta idxIabc]);

% Equation for TorqueM
% Depends on gov control, so leave row 17 of Jacobian empty for now!

idx = idx(1:idxCount);
jdx = jdx(1:idxCount);

end