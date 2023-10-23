function matP = matPolytope(matZ)
% matPolytope - Converts a matrix zonotope into a matrix polytope 
%    representation
%
% Syntax:
%    matP = matPolytope(matZ)
%
% Inputs:
%    matZ - matZonotope object
%
% Outputs:
%    matP - matPolytope object
%
% Example: 
%
% Other m-files required: vertices, polytope
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff
% Written:       22-June-2010
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% obtain vertices
V=vertices(zonotope(matZ));

%obtain vertices
for i=1:length(V(1,:))
    matrixVertex{i}=vec2mat(V(:,i));
end

%convert polytope to matrix polytope
matP=matPolytope(matrixVertex);

% ------------------------------ END OF CODE ------------------------------
