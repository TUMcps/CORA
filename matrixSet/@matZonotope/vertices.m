function matV = vertices(matZ)
% vertices - computes the vertices of a matrix zonotope
%
% Syntax:  
%    matV = vertices(matZ)
%
% Inputs:
%    matZ - matrix zonotope
%
% Outputs:
%    matV - cell array of matrix vertices
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      24-June-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%conversion to an intervalhull
Z = zonotope(matZ);

%compute vertices
V = vertices(Z,'polytope');
V = unique(V', 'rows')'; %eliminate vectors that occur multiple times

%convert vertices to matrix vertices
matV=cell(length(V(1,:)),1);
for i=1:length(V(1,:))
    matV{i}=vec2mat(V(:,i));
end

%------------- END OF CODE --------------