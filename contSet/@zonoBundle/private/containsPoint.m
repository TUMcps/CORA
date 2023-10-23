function res = containsPoint(Z,p)
% containsPoint - checks if p inside zonoBundle
%
% Syntax:
%    res = containsPoint(Z,p)
%
% Inputs:
%    Z - zonoBundle object
%    p - point(s) specified as a column vector (array)
%
% Outputs:
%    res - boolean (array) whether p in Z
%
% Example: 
%    Z1 = zonotope([1;1], [1 1; -1 1]);
%    Z2 = zonotope([-1;1], [1 0; 0 1]);
%    zB = zonoBundle({Z1,Z2});
%    p = [-0.5;1];
%    contains(zB,p)
% 
%    figure; hold on;
%    plot(zB); scatter(p(1),p(2),8,'r','filled');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       18-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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

% ------------------------------ END OF CODE ------------------------------
