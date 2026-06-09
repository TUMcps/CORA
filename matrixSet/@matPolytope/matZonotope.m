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
%                13-April-2026 (LS, fix shape of matrix zonotope)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get (matrix) vertices
matV = matP.V;

%convert matrix vertices to vertices
V = reshape(matV,[],size(matV,3));

%convert to zonotope
Z = compact(zonotope.enclosePoints(V),'zeros');

% reshape center and generator matrix
C = reshape(Z.c,[size(matP),1]);
G = reshape(Z.G,[size(matP),size(Z.G,2)]);

%convert to matrix zonotope
matZ = matZonotope(C,G);

% ------------------------------ END OF CODE ------------------------------
