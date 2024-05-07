function intMat = intervalMatrix(matZ)
% intervalMatrix - computes an enclosing interval matrix of a matrix
%    zonotope
%
% Syntax:
%    intMat = intervalMatrix(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    intMat - intervalMatrix object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       21-June-2010 
% Last update:   06-October-2010
%                26-August-2011
%                03-April-2023 (MW, remove setting property)
%                25-April-2024 (TL, much faster implementation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% center matrix
C = matZ.C;

% compute delta
D = sum(abs(matZ.G),3);

% instantiate interval matrix
intMat = intervalMatrix(C, D);

end

% ------------------------------ END OF CODE ------------------------------
