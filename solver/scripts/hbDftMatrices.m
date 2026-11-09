function [GAMMA,invGAMMA] = hbDftMatrices(N)
%HBDFTMATRICES  Real-valued DFT pair for the harmonic-balance Jacobian.
%
%   [GAMMA,invGAMMA] = hbDftMatrices(N) builds the transform pair used by
%   every nonlinear device stamp:
%
%       X = GAMMA*x        time samples -> frequency coefficients
%       x = invGAMMA*X     frequency coefficients -> time samples
%
%   GAMMA is [2*ceil(N/2) x N] and invGAMMA is [N x 2*ceil(N/2)].
%
%   Two conventions are baked in:
%
%   1. The signals are REAL, so the negative-frequency coefficients are the
%      complex conjugates of the positive ones and carry no new
%      information. Only DC and the positive frequencies are kept. The
%      factor of 2 on invGAMMA restores the conjugate partner's
%      contribution; DC has no partner, so its two columns are halved back.
%
%   2. Complex coefficients are stored as interleaved real and imaginary
%      parts, [Xre(1); Xim(1); Xre(2); Xim(2); ...], which is the same
%      (real,imag) interleaving the rest of the solver uses. That is what
%      makes GAMMA*diag(d)*invGAMMA a real matrix and lets it be stamped
%      directly.
%
%   Row 2 of GAMMA - the imaginary part of the DC coefficient - is
%   identically zero, because a real signal has a real DC term. That is
%   what makes every nonlinear block's imaginary-DC row structurally zero,
%   and why HB.buildStampMaps has to supply those rows explicitly. See
%   HB.stampLinear_Y_I.
%
%   NOTE: dftmtx(N) would do the first line, but it lives in the Signal
%   Processing Toolbox, which is not installed everywhere. MathWorks
%   documents dftmtx(n) as fft(eye(n)), so this is the same matrix with no
%   toolbox dependency.

numFreqVars = 2*ceil(N/2);      % DC + positive frequencies, each (re,im)
idxDC = floor(N/2)+1;           % DC row after the fftshift

% ----------------------- Time to frequency ------------------------------
GAMMA_complex = fftshift(fft(eye(N)),1);

GAMMA = zeros(numFreqVars,N);
GAMMA(1:2:end,:) = real(GAMMA_complex(idxDC:end,:));
GAMMA(2:2:end,:) = imag(GAMMA_complex(idxDC:end,:));

% ----------------------- Frequency to time ------------------------------
invGAMMA_complex = inv(GAMMA_complex);                          %#ok<MINV>

% Re{ (Gre + j*Gim)*(Xre + j*Xim) } = Gre*Xre - Gim*Xim
invGAMMA = zeros(N,numFreqVars);
invGAMMA(:,1:2:end) =  2*real(invGAMMA_complex(:,idxDC:end));
invGAMMA(:,2:2:end) = -2*imag(invGAMMA_complex(:,idxDC:end));
% DC has no conjugate partner
invGAMMA(:,1:2) = invGAMMA(:,1:2) / 2;

end
