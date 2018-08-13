function eZ = expmVertex(matZ)
% expmVertex - computes the exponential matrix for the vertices of the
% matrix zonotope
%
% Syntax:  
%    eZ = expmVertex(matZ)
%
% Inputs:
%    matZ - matrix zonotope
%
% Outputs:
%    eZ - matrix zonotope exponential
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---end+i

% Author:       Matthias Althoff
% Written:      19-July-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%compute vertices
V = dominantVertices(matZ, 100);

%compute exponential matrices
for i=1:length(V)
    eZ{i} = expm(V{i});
end


%------------- END OF CODE --------------