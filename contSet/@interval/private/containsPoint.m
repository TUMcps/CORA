function res = containsPoint(Int,p)
% containsPoint - determines if the point p is inside the interval Int
%
% Syntax:  
%    result = containsPoint(Int,p)
%
% Inputs:
%    Int - interval object
%    p - point(s) specified as a column vector (array)
%
% Outputs:
%    result - boolean (array) whether point is inside the interval or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

numPoints = size(p,2);
res = false(numPoints,1);
infi = infimum(Int);
supr = supremum(Int);

for iPoint = 1:numPoints
    p_curr = p(:,iPoint);
    res(iPoint) = all(infi <= p_curr) && all(supr >= p_curr);
end

%------------- END OF CODE --------------
