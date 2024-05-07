function matZ = matZonotope(matP)
% matZonotope - computes an enclosing matrix zonotope of a matrix polytope
%
% Syntax:
%    matZ = matZonotope(matP)
%
% Inputs:
%    matP - matPolytope object
%
% Outputs:
%    matZ - matZonotope object
%
% Example: 
%    - 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Authors:       Matthias Althoff
% Written:       22-July-2010 
% Last update:   02-May-2024 (TL, new structure of vertices)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get (matrix) vertices
matV = matP.V;

%convert matrix vertices to vertices
V = reshape(matV,[],size(matV,3));

%convert to zonotope
Z = zonotope.enclosePoints(V);

%convert to matrix zonotope
matZ = matZonotope(Z);

% ------------------------------ END OF CODE ------------------------------
