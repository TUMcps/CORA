function res = isequal(Z1,Z2)
% isequal - checks if two zonoBundles are equal
%
% Syntax:  
%    res = isequal(Z1,Z2)
%
% Inputs:
%    Z1 - zonoBundle object
%    Z2 - zonoBundle object
%
% Outputs:
%    res - boolean whether Z1 and Z2 are equal
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
if Z1.parallelSets ~= Z2.parallelSets
    return
end

for z=1:Z1.parallelSets
    if ~isequal(Z1.Z{z},Z2.Z{z})
        return
    end
end

res = true;

%------------- END OF CODE --------------