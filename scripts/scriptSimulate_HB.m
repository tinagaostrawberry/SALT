%% Newton-Raphson on the HB harmonic-balance system.
%
% The same loop as SCRIPTSIMULATE_SALT - refill the triplet values, rebuild
% the sparse Jacobian, take one step of
% J(x[k])*x[k+1] = J(x[k])*x[k] - F(x[k]) - with one addition that is
% specific to harmonic balance: the fundamental w0 is an UNKNOWN, so the
% sample grid and the DFT matrices are rebuilt from the new w0 at the end of
% every step. That rebuild is also what makes the Park transforms free of
% any d/dw0 term, because it keeps w0*tSamp(n) = 2*pi*n/N whatever w0 is.
%
% Reads   hb, totIters_NR, relErrTol_NR, absErrTol_NR,
%         include_clampLimiting, clampLimit
% Writes  hb, converge_HB, idx_HB, time_HB

idx_HB = 0;
converge_HB = false;
tic
while (idx_HB < totIters_NR) && (~converge_HB)

    % Stamp
    hb = hb.stampLinear_Y_I;
    hb = hb.stampNonlinear_Y_I;
    hb = hb.stampPhase_Y_I;
    hb = hb.assemble_Y_I;

    % Solve
    x_v_prevGuess = hb.x_v;
    hb.x_v = hb.Y \ hb.I;

    % Check convergence BEFORE any limiting, so the heuristic cannot make
    % the step look converged when it is not
    converge_HB = isErrSmall(hb.x_v,x_v_prevGuess,relErrTol_NR,absErrTol_NR);

    % Absolute limiting - voltages are PU and the frequency should not move
    % far, so no single step should exceed clampLimit.
    %
    % SCALED BY N, which is not a fudge. clampLimit is a limit on a PHYSICAL
    % quantity, and x_v holds UNNORMALIZED DFT coefficients: a DC physical
    % value c appears as c*N (and a phasor as c*N/2). Clamping the
    % coefficients at clampLimit would therefore hold HB to 1/N of the step
    % SALT is allowed, which is exactly what it did before this line was
    % scaled - HB needed two to three times SALT's iterations on every
    % contingency, in lockstep with N = 3, and the hardest contingency ran
    % out of iterations altogether.
    if include_clampLimiting
        hb.x_v = stepLimiting(hb.x_v,x_v_prevGuess,[0 hb.N*clampLimit]);
    end

    % Every tone is k*w0, so the sample grid and both DFT matrices follow
    % the frequency
    [hb.f0,hb.tSamp,hb.fSamp,hb.N,hb.GAMMA,hb.invGAMMA,hb.fGamma] = ...
        hb.updateFreqSamples(hb.x_v(hb.getRowIdx_w0));
    [hb.tSim,hb.deltaT] = hb.updateTimeSimSamples;

    idx_HB = idx_HB + 1;
end
time_HB = toc;

if ~converge_HB
    error(['HB did not converge after ' num2str(idx_HB) ' iteration(s)'])
end
disp(['HB converged after ' num2str(idx_HB) ' iteration(s).'])
disp(['HB solve time: ' num2str(time_HB) ' seconds.'])
