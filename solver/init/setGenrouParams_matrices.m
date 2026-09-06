function [RMatrix,LMatrix,LspMatrix] = setGenrouParams_matrices( ...
    genrou_Ra,genrou_Rfd,genrou_R1d,genrou_R1q,genrou_R2q, ...
    genrou_Lad,genrou_Laq,genrou_Ll,genrou_L0, ...
    genrou_Lffd,genrou_Lf1d,genrou_L11d,genrou_L11q,genrou_L22q)
%SETGENROUPARAMS_MATRICES  GENROU resistance and inductance matrices.
%
%   Every input is a column vector with one entry per machine. Returns three
%   [7 x 7 x numGenrou] arrays, one page per machine, in the winding order
%
%       [ d  q  0 | fd  1d | 1q  2q ]
%
%   RMatrix    winding resistances (diagonal). The rotor entries are
%              negative because of the equivalent-circuit orientation.
%   LMatrix    self and mutual inductances.
%   LspMatrix  the flux terms multiplying rotor speed in the stator
%              equation: the d and q rows exchanged and signed.

col = @(v) v(:);

% Resistance: diagonal
rowR = 1:7;
colR = 1:7;
datR = [repmat(col(genrou_Ra),1,3), ...
        -col(genrou_Rfd), -col(genrou_R1d), ...
        -col(genrou_R1q), -col(genrou_R2q)];

% Inductance
%          phid   phiq  phi0 | phifd phi1d | phi1q  phi2q
rowL = [1 1 1, 2 2 2, 3, 4 4 4, 5 5 5, 6 6 6, 7 7 7];
colL = [1 4 5, 2 6 7, 3, 1 4 5, 1 4 5, 2 6 7, 2 6 7];
datL = [-(col(genrou_Lad)+col(genrou_Ll)),  col(genrou_Lad),  col(genrou_Lad), ...
        -(col(genrou_Laq)+col(genrou_Ll)),  col(genrou_Laq),  col(genrou_Laq), ...
        -col(genrou_L0), ...
        -col(genrou_Lad), col(genrou_Lffd), col(genrou_Lf1d), ...
        -col(genrou_Lad), col(genrou_Lf1d), col(genrou_L11d), ...
        -col(genrou_Laq), col(genrou_L11q), col(genrou_Laq), ...
        -col(genrou_Laq), col(genrou_Laq),  col(genrou_L22q)];

% Speed-dependent flux terms
rowS = [1 1 1, 2 2 2];
colS = [2 6 7, 1 4 5];
datS = [ (col(genrou_Laq)+col(genrou_Ll)), -col(genrou_Laq), -col(genrou_Laq), ...
        -(col(genrou_Lad)+col(genrou_Ll)),  col(genrou_Lad),  col(genrou_Lad)];

RMatrix   = scatterPages(rowR,colR,datR);
LMatrix   = scatterPages(rowL,colL,datL);
LspMatrix = scatterPages(rowS,colS,datS);
end

function M = scatterPages(row,col,dat)
% dat is [numGenrou x numel(row)]: one row per machine, one column per
% nonzero entry of the pattern.
n = size(dat,1);
M = zeros(7,7,n);
linIdx = sub2ind([7,7],row,col);
for k = 1:n
    page = zeros(7);
    page(linIdx) = dat(k,:);
    M(:,:,k) = page;
end
end
