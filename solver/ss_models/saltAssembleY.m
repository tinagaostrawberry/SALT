function Y = saltAssembleY(row,col,val,n)
%SALTASSEMBLEY  Build the sparse Jacobian from its (value,row,col) triplets.
%
%   Duplicate (row,col) pairs are summed, which is exactly the accumulate
%   semantics the network stamps rely on. Structural zeros are dropped.

Y = sparse(row,col,val,n,n);
end
