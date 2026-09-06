function x = limitAbsVal(x,minMax)
% Adjust to minimum absolute limit
x(abs(x)<minMax(1)) = sign(x(abs(x)<minMax(1)))*minMax(1);
% Adjust to maximum absolute limit
x(abs(x)>minMax(2)) = sign(x(abs(x)>minMax(2)))*minMax(2);
end