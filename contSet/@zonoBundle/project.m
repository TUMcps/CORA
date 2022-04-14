function Zbundle = project(Zbundle,dims)
% project - Returns a zonotope bundle which is projected onto the specified
% dimensions
%
% Syntax:  
%    Zbundle = project(Zbundle,dim)
%
% Inputs:
%    Zbundle - zonotope bundle
%    dims - vector of two dimensions onto which one should project
%
% Outputs:
%    Zbundle - zonotope bundle
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      04-February-2011
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%project for each zonotope
for i=1:Zbundle.parallelSets
    Zbundle.Z{i}=project(Zbundle.Z{i},dims);
end

%------------- END OF CODE --------------