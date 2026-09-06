function x_v_updatedGuess_limited = stepLimiting(x_v_updatedGuess,x_v_prevGuess, ...
    LIMITS)

dx_og = x_v_updatedGuess - x_v_prevGuess;
dx = limitAbsVal(dx_og,LIMITS);

% Apply limited step
x_v_updatedGuess_limited = x_v_prevGuess + dx;

end