function [idx,jdx, ...
    idx_dFx2_dwr,jdx_dFx2_dwr, ...
    idx_dFx3_dwr,jdx_dFx3_dwr, ...
    idx_dFTm_dPmech,jdx_dFTm_dPmech, ...
    idx_dFTm_dwr,jdx_dFTm_dwr, ...
    idx_dFTm_dTm,jdx_dFTm_dTm] = calcTgov1dIndices


% Indices
[~,~,~,idxWr,~,~,idxTm] = getGenrouIndicesTD;

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

% State 1 (backpropgation for x2)
% x1 + T3*dx1/dt - x2 - T2*dx2/dt = 0
append_block(1,[1 2]);

% State 2 (backpropagation for pref)
% pref - (wr/ws-1) - R*x2 - R*T1*dx2/dt = 0
% =>
% pref - ((wr+wrPrev)/2)/ws + 1 - R*(x2+x2Prev)/2 - R*T1*(x2-x2Prev)/T = 0
append_block(2,2);
%
[idx_dFx2_dwr,jdx_dFx2_dwr] = ndgrid(2,idxWr);

% Pmech/state 3
% x1 - (wr/ws-1)*Dt - pmech = 0
append_block(3,[1 3]);
%
[idx_dFx3_dwr,jdx_dFx3_dwr] = ndgrid(3,idxWr);

idx = idx(1:idxCount);
jdx = jdx(1:idxCount);

% Check electrical and mechanical angles are equal
% theta = (p/2)*thetaM which means
% wr = (p/2)*wM, since w = dtheta/dt. If p = 2, then wr = (2/2)*wM = wM.
% When wr = wM, we have:
% pmech - (wr/ws)*Tm = 0
[idx_dFTm_dPmech,jdx_dFTm_dPmech] = ndgrid(idxTm,3);
[idx_dFTm_dwr,jdx_dFTm_dwr] = ndgrid(idxTm,idxWr);
[idx_dFTm_dTm,jdx_dFTm_dTm] = ndgrid(idxTm,idxTm);


end