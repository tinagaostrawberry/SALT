function [idx,jdx, ...
    idx_dFefd_dx1,jdx_dFefd_dx1, ...
    idx_dFefd_defd,jdx_dFefd_defd, ...
    idx_dFefd_dwr,jdx_dFefd_dwr, ...
    idx_dFexcVt_dvt,jdx_dFexcVt_dvt, ...
    idx_dFexcVt_dvabc,jdx_dFexcVt_dvabc] = calcExcDc1Indices


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

% State 1 (backpropgation for x4)
% -Ke*x1(t) - Te*dx1(t)/dt + x4(t) = 0
append_block(1,[1 4]);

% State 2 (backpropagation for x2)
% x2(t) + Tr*dx2(t)/dt - vt(t) = 0
append_block(2,[2 6]);

% State 3 (backpropagation for vref)
% x2(t) + Tc*dx2(t)/dt +
% x3(t) + Tb*dx3(t)/dt +
% x5(t) + Tc*dx5(t)/dt +
% - vref(t) - Tc*dvref(t)/dt = 0
%            |______________|
%                   0
append_block(3,[2 3 5]);

% State 4 (backpropagation for x3)
% -Ka*x3(t) + x4(t) + Ta*dx4(t)/dt = 0
append_block(4,[3 4]);

% State 5 (backpropagation for x5)
% -Kf1*dx1(t)/dt + x5(t) + Tf1*dx5(t)/dt = 0
append_block(5,[1 5]);

idx = idx(1:idxCount);
jdx = jdx(1:idxCount);

% Efd normalized by speed
% excDc1_efd - x1*wr_pu = 0
% excDc1_efd = genToExc*genrou_efd
% =>
% genToExc*genrou_efd - x1*wr_pu = 0
[~,~,idxEfd,idxWr] = getGenrouIndicesTD;
[idx_dFefd_dx1,jdx_dFefd_dx1] = ndgrid(idxEfd,1);
[idx_dFefd_defd,jdx_dFefd_defd] = ndgrid(idxEfd,idxEfd);
[idx_dFefd_dwr,jdx_dFefd_dwr] = ndgrid(idxEfd,idxWr);

% Vt
% -Vt + sqrt(2/3)*(va^2 + vb^2 + vc^2)^(1/2) = 0 s.t. va,vb,vc are gen PU
[idx_dFexcVt_dvt,jdx_dFexcVt_dvt] = ndgrid(6,6);
[idx_dFexcVt_dvabc,jdx_dFexcVt_dvabc] = ndgrid(6,[0 0 0]);
idx_dFexcVt_dvabc = idx_dFexcVt_dvabc(:);
jdx_dFexcVt_dvabc = jdx_dFexcVt_dvabc(:);

end