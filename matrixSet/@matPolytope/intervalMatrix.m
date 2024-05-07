function intMat = intervalMatrix(matP)
% intervalMatrix - computes an enclosing interval matrix of a matrix
%    polytope
%
% Syntax:
%    intMat = intervalMatrix(matP)
%
% Inputs:
%    matP - matPolytope object
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
% Last update:   02-May-2024 (TL, simplified due to new structure of V)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain vertices
V = matP.V;

% compute min/max across all vertices
minMat= min(V,[],3);
maxMat= max(V,[],3);

% instantiate interval matrix
C=0.5*(minMat+maxMat);
D=0.5*(maxMat-minMat);
intMat=intervalMatrix(C,D);

% ------------------------------ END OF CODE ------------------------------
