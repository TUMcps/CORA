function matV = dominantVertices(matZ, maxNumber)
% dominantVertices - computes the dominant vertices of a matrix zonotope
%
% Syntax:  
%    matV = dominantVertices(matZ, maxNumber)
%
% Inputs:
%    matZ - matrix zonotope
%    maxNumber - maximum number of dominant vertices
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
% Written:      19-July-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%conversion to a zonotope
Z = zonotope(matZ);

%underapproximate zonotope
V = underapproximate(Z);

%convert vertices to matrix vertices
for i=1:length(V(1,:))
    matV{i}=vec2mat(V(:,i));
end

%------------- END OF CODE --------------