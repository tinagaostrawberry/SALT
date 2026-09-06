function Y_ReIm = stampReIm(Y)
% Implements:
% [real(Y) -imag(Y); imag(Y) real(Y)];


Y_ReIm = zeros(2*size(Y));

Y_ReIm(1:2:end,1:2:end) = real(Y);
Y_ReIm(1:2:end,2:2:end) = -imag(Y);
Y_ReIm(2:2:end,1:2:end) = imag(Y);
Y_ReIm(2:2:end,2:2:end) = real(Y);


end