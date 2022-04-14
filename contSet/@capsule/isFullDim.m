function res = isFullDim(obj)
% isFullDim - check if a capsule is full-dimensional
%
% Syntax:  
%    res = isFullDim(obj)
%
% Inputs:
%    obj - capsule object
%
% Outputs:
%    res - true if capsule is full-dimensional, otherwise false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/isFullDim

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      02-January-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = true;
if obj.r == 0
    % just a line in this case
    res = false;
end

%------------- END OF CODE --------------