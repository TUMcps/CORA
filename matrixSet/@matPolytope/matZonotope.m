function matZ = matZonotope(matP)
% matZonotope - computes an enclosing matrix zonotope of a matrix
% polytope
%
% Syntax:  
%    matZ = matZonotope(matP)
%
% Inputs:
%    matP - matrix polytope
%
% Outputs:
%    matZ - matrix zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      22-July-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%convert matrix polytope to matrix vertices
matV = vertices(matP);

%convert matrix vertices to vertices
V=[];
for i=1:length(matV)
    V(:,end+1)=mat2vec(matV{i});
end

%convert to zonotope
Z = zonotope(vertices(V));

%convert to matrix zonotope
matZ = matZonotope(Z);

%------------- END OF CODE --------------