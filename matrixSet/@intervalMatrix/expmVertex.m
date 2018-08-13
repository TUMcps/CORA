function eI = expmVertex(intMat)
% expmVertex - computes the exponential matrix for the vertices of the
% interval matrix
%
% Syntax:  
%    eI = expmVertex(intMat)
%
% Inputs:
%    intMat - interval matrix
%
% Outputs:
%    eI - interval matrix exponential
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---end+i

% Author:       Matthias Althoff
% Written:      02-July-2010 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%compute vertices
V = dominantVertices(intMat, 100);

%compute exponential matrices
for i=1:length(V)
    eI{i} = expm(V{i});
end


%------------- END OF CODE --------------