function x = getTimeDomainSignal_iDFT(tSamp,funGetIndices, ...
    x_v,frequencies,K_tot)

numT = numel(tSamp);

% Solve at DC
rowIndices = funGetIndices(0);
if any(frequencies==0)
    % Get DC coefficient
    X = x_v(rowIndices);
    assert(all(abs(X(2:2:end))<1e-9))
    x = X(1:2:end) .* ones(1,numT);
    % Remove DC component going forward
    frequencies(frequencies==0) = [];
    K_tot = K_tot-1;
    assert(numel(frequencies)==K_tot)
else
    x = zeros(numel(rowIndices)/2,numT);
end

% Solve at other frequencies
for idxFreq = 1:K_tot
    wk = 2*pi*frequencies(idxFreq);

    rowIndices = funGetIndices(idxFreq);
    X = x_v(rowIndices);
    % Assume real, so also include complex conjugate, which
    % results in 2x real part
    % x = x + X.*exp(1j*w*tSamp) + (X').*exp(-1j*w*tSamp);
    x = x + 2*real( complex(X(1:2:end),X(2:2:end)).*exp(1j*wk*tSamp.') );
end

% Normalize by number data points, as convention of iDFT
x = x/numT;

end