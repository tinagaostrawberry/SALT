function X_ReIm = hbDft(GAMMA,K,x)
%HBDFT  Harmonic coefficients of a set of time-domain signals.
%
%   X_ReIm = hbDft(GAMMA,K,x) takes x, [nVar x N] time samples, and returns
%   [2*nVar x (K+1)]: one column per harmonic 0..K, with the rows carrying
%   each variable's real and imaginary part interleaved.
%
%   Flattening the result column-major, X_ReIm(:), therefore runs
%   variable-fastest within each harmonic - which is exactly the row order
%   HB uses inside a device block (k*nEqns + localIndex). That is what lets
%   a kernel hand back a flat residual and have it land correctly.
%
%   This is the forward direction of the pair HBDFTMATRICES builds; HBIDFT
%   is the inverse. Note the two are NOT normalised the same way: x_v holds
%   unnormalized coefficients, so GAMMA carries no 1/N and HBIDFT does.

X = (GAMMA*(x.')).';                    % nVar x 2*(K+1)

X_ReIm = zeros(size(x,1)*2,K+1);
X_ReIm(1:2:end,:) = X(:,1:2:end);
X_ReIm(2:2:end,:) = X(:,2:2:end);

end
