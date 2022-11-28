function res = isequal(zB1,zB2)
% isequal - checks if two zonotope bundles are equal
%
% Syntax:  
%    res = isequal(zB1,zB2)
%
% Inputs:
%    zB1 - zonoBundle object
%    zB2 - zonoBundle object
%
% Outputs:
%    res - true/false
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = false;
if zB1.parallelSets ~= zB2.parallelSets
    return
end

for z=1:zB1.parallelSets
    if ~isequal(zB1.Z{z},zB2.Z{z})
        return
    end
end

res = true;

%------------- END OF CODE --------------