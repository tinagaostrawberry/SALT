%% Newton-Raphson on the SALT steady-state system.
%
% Each iteration refills the triplet values (no element-wise writes into a
% dense matrix), rebuilds the sparse Jacobian and takes one Newton step in
% the form J(x[k])*x[k+1] = J(x[k])*x[k] - F(x[k]).

idx_SS = 0;
converge_SS = false;
tic
while (idx_SS < totIters_NR) && (~converge_SS)

    % Stamp
    salt = salt.stampLinear_Y_I;
    salt = salt.stampNonlinear_Y_I;
    salt.stampPhase_Y_I;
    salt = salt.assemble_Y_I;

    % Solve
    x_v_prevGuess = salt.x_v;
    salt.x_v = salt.Y \ salt.I;

    % Check convergence BEFORE any limiting, so the heuristic cannot make
    % the step look converged when it is not
    converge_SS = isErrSmall(salt.x_v,x_v_prevGuess,relErrTol_NR,absErrTol_NR);
    [~,maxIdx] = max(abs(salt.x_v-x_v_prevGuess));
    disp(['Max err (@ idx=' num2str(maxIdx) '): ' ...
        num2str(salt.x_v(maxIdx)-x_v_prevGuess(maxIdx))])

    % Absolute limiting - voltages are PU and the frequency should not
    % move far, so no single step should exceed 1
    if include_clampLimiting
        salt.x_v = stepLimiting(salt.x_v,x_v_prevGuess,[0 clampLimit]);
    end

    % The period (and therefore the sample grid) follows the frequency
    salt = salt.updateTimeSimSamples;

    idx_SS = idx_SS + 1;
end
time_SS = toc;

if ~converge_SS
    error(['SALT did not converge after ' num2str(idx_SS) ' iteration(s)'])
end
disp(['SALT converged after ' num2str(idx_SS) ' iteration(s).'])
disp(['SALT solve time: ' num2str(time_SS) ' seconds.'])
