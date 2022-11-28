function res = isFullDim(C)
% isFullDim - checks if the dimension of the affine hull of a capsule is
%    equal to the dimension of its ambient space
%
% Syntax:  
%    res = isFullDim(C)
%
% Inputs:
%    C - capsule object
%
% Outputs:
%    res - true/false
%
% Example:
%    C = capsule([2;1],[1;1],1);
%    isFullDim(C)
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
if C.r == 0
    % just a point/line in this case
    res = false;
end

%------------- END OF CODE --------------