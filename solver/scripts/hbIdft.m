function x = hbIdft(tSamp,funGetIndices,x_v,frequencies,K_tot)
%HBIDFT  Reconstruct time samples from harmonic coefficients in x_v.
%
%   x = hbIdft(tSamp,funGetIndices,x_v,frequencies,K_tot)
%
%   funGetIndices(k) must return the rows of x_v holding one variable set
%   at harmonic k, as interleaved (real,imag) pairs. The returned x is
%   [nVar x numel(tSamp)].
%
%   The signals are real, so each positive-frequency coefficient is added
%   together with its conjugate partner - hence the 2*real(...). DC has no
%   partner and is added once.
%
%   The 1/numel(tSamp) at the end is the DFT convention this solver uses:
%   x_v holds UNNORMALIZED coefficients, i.e.
%       x(n) = (1/N) * sum_k X(k)*exp(1j*2*pi*k*n/N).
%   Every scaling in HB follows from that one choice - see
%   HB.systemValInitializations.

numT = numel(tSamp);

% DC first, if it is in the set
rowIndices = funGetIndices(0);
if any(frequencies==0)
    X = x_v(rowIndices);
    % A real signal has a real DC coefficient. HB enforces this
    % structurally, so a nonzero value here means the enforcement broke.
    assert(all(abs(X(2:2:end))<1e-9), ...
        'hbIdft: DC coefficient has a nonzero imaginary part')
    x = X(1:2:end) .* ones(1,numT);
    frequencies(frequencies==0) = [];
    K_tot = K_tot-1;
    assert(numel(frequencies)==K_tot)
else
    x = zeros(numel(rowIndices)/2,numT);
end

% Then every retained harmonic
for idxFreq = 1:K_tot
    wk = 2*pi*frequencies(idxFreq);
    X = x_v(funGetIndices(idxFreq));
    x = x + 2*real( complex(X(1:2:end),X(2:2:end)).*exp(1j*wk*tSamp.') );
end

x = x/numT;

end
