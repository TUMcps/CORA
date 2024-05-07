function matZ = concatenate(matZ1,matZ2)
% concatenate - concatenates the center and all generators of the second
%    matrix zonotope to the first one
%
% Syntax:
%    matZ = concatenate(matZ1,matZ2)
%
% Inputs:
%    matZ1 - matZonotope object
%    matZ2 - matZonotope object
%
% Outputs:
%    matZ - matrix zonotope object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       06-September-2013
% Last update:   25-April-2024 (TL, new implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input
if ~all(dim(matZ1) == dim(matZ2)) ...
    && ~(isempty(matZ1) || isempty(matZ2))
    throw(CORAerror('CORA:dimensionMismatch',matZ1,matZ2))
end

% read properties
C1 = matZ1.C; G1 = matZ1.G;
C2 = matZ2.C; G2 = matZ2.G;

% init new matrix zonotope
C = C1;
G = cat(3,G1,C2,G2);
matZ = matZonotope(C,G);

end

% ------------------------------ END OF CODE ------------------------------
