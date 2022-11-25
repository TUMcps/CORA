function res = containsPoint(Z,p)
% containsPoint - checks if p inside zonoBundle
%
% Syntax:  
%    res = containsPoint(Z)
%
% Inputs:
%    Z - zonoBundle object
%    p - point(s) specified as a column vector (array)
%
% Outputs:
%    res - boolean (array) whether p in Z
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
% Written:      18-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

numPoints = size(p,2);
res = false(numPoints,1);

for iPoint = 1:numPoints
    p_curr = p(:,iPoint);
    
    for z=1:Z.parallelSets
        if ~containsPoint(Z.Z{z},p_curr)
            res(iPoint) = false;
            break
        end
    end
end

%------------- END OF CODE --------------