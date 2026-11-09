function [idxReIm,idxFreq] = hbSelectDftHarmonics(fGamma,fInterest)
%HBSELECTDFTHARMONICS  Rows/columns of the DFT pair for the kept harmonics.
%
%   [idxReIm,idxFreq] = hbSelectDftHarmonics(fGamma,fInterest)
%
%   fGamma     the non-negative frequencies HBDFTMATRICES built, in order
%   fInterest  the harmonics the solve retains, 0 .. K times f0
%
%   idxFreq   index of each retained harmonic within fGamma
%   idxReIm   the same, expanded to the interleaved (real,imag) row/column
%             pairs that GAMMA and invGAMMA actually carry
%
%   With N = 2*K+1 the two sets coincide and this is the identity. It
%   exists so that oversampling the time grid - which is how a harmonic
%   balance normally controls aliasing in the Jacobian - needs no change
%   anywhere else. See HBTIMESAMPLES.

numReImEqns = 2;

[found,idxFreq] = ismembertol(fInterest,fGamma);
assert(all(found), ...
    'hbSelectDftHarmonics: a retained harmonic is not on the DFT grid')

idxReIm = numReImEqns*reshape(idxFreq-1,1,[]) + (1:numReImEqns).';
idxReIm = idxReIm(:);

end
