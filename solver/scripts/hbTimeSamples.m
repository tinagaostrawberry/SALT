function [tSamp,fSamp] = hbTimeSamples(K,f0)
%HBTIMESAMPLES  Sample grid for a single-tone harmonic-balance solve.
%
%   [tSamp,fSamp] = hbTimeSamples(K,f0) returns N = 2*K+1 time samples
%   spanning exactly one period of f0, and the matching two-sided frequency
%   grid (-K:K)*f0.
%
%   Every tone is k*f0, so the harmonics are commensurate by construction
%   and an ordinary DFT suffices. Incommensurate tones (an interharmonic
%   from a converter, say) would need an almost-periodic transform instead
%   and are not representable in this formulation at all.
%
%   NOTE ON OVERSAMPLING. N = 2*K+1 is the MINIMUM grid that resolves K
%   harmonics, so there is no anti-aliasing margin: in the Jacobian the
%   DFT sandwich folds everything a nonlinearity generates above K back
%   onto the retained harmonics. Standard harmonic-balance practice
%   oversamples. The machinery for that already exists - HBDFTMATRICES
%   builds the transform on whatever grid it is handed, and
%   HB.updateFreqSamples then selects only the harmonics of interest - so
%   raising N here is the single change needed. It is left at the minimum
%   to match the legacy solver exactly.

tolErr = 1e-12;

% Frequency grid, two-sided
fSamp = [(-K:1:-1) 0 (1:K)].'*f0;

% Time grid: one period of f0, N uniform samples
freqBinSize = f0;
numFreqs = numel(fSamp);
deltaT = 1/(freqBinSize*numFreqs);
tSamp = (0:deltaT:(deltaT*(numFreqs-1))).';

% The three properties the harmonic-balance derivation relies on
f_set = (0:K).'*f0;
% 1. every retained harmonic is recoverable from a DFT on this grid
assert(all(ismembertol(f_set,hbFreqFromTime(tSamp),tolErr)), ...
    'HB sample grid does not resolve every retained harmonic')
% 2. the grid spans a whole number of periods of every retained harmonic,
%    so the DFT sees no leakage
multPeriod = (tSamp(end)+deltaT)./(1./f_set);
assert(all(abs(multPeriod-round(multPeriod))<tolErr), ...
    'HB sample grid is not periodic in every retained harmonic')
% 3. no aliasing of the retained harmonics themselves
assert(deltaT<1/(2*max(f_set)), ...
    'HB sample grid aliases the highest retained harmonic')

end
