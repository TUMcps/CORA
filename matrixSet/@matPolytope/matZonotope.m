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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%convert matrix polytope to matrix vertices
matV = vertices(matP);

%convert matrix vertices to vertices
V=[];
for i=1:length(matV)
    V(:,end+1)=mat2vec(matV{i});
end

%convert to zonotope
Z = zonotope.enclosePoints(V);

%convert to matrix zonotope
matZ = matZonotope(Z);

% ------------------------------ END OF CODE ------------------------------
