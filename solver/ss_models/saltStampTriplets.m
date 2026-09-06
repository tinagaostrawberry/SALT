function [row,col] = saltStampTriplets(rowIdx,colIdx)
%SALTSTAMPTRIPLETS  Expand one dense Jacobian block into (row,col) triplets.
%
%   [row,col] = saltStampTriplets(rowIdx,colIdx) returns the global row and
%   column index of every entry of the dense block A that a scalar solver
%   would write as  Y(rowIdx,colIdx) = A.
%
%   The expansion is column-major, so the returned lists line up entry for
%   entry with A(:). That is what lets the vectorised kernels hand back a
%   flat value array and have it land in the right places.

[R,C] = ndgrid(rowIdx(:),colIdx(:));
row = R(:);
col = C(:);
end
