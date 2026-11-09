function fSamp = hbFreqFromTime(tSamp)
%HBFREQFROMTIME  Two-sided DFT frequency grid of a uniform time grid.
%
%   fSamp = hbFreqFromTime(tSamp) returns the frequencies, in the same
%   fftshift-ed (-N/2 .. N/2-1) order HBDFTMATRICES uses, that an N-point
%   DFT on the uniform grid tSamp resolves.
%
%   Reference: https://gubner.ece.wisc.edu/notes/FrequencyAnalysisLong.pdf

N = numel(tSamp);
deltaT = uniquetol(diff(tSamp),1e-9);
assert(isscalar(deltaT),'hbFreqFromTime needs a uniform time grid')

fs = 1/deltaT;
kvec = (0:N-1) - floor(N/2);
fSamp = ((kvec/N)*fs).';

end
