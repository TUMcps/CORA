function res = containsPoint(I,p)
% containsPoint - determines if an interval contains a point cloud p, where
%    the result is evaluated for each point
%
% Syntax:  
%    result = containsPoint(I,p)
%
% Inputs:
%    I - interval object
%    p - point(s) specified as a column vector (array)
%
% Outputs:
%    result - true/false (array)
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

% check input arguments
inputArgsCheck({{I,'att','interval'};
                {p,'att','numeric','column'}});

numPoints = size(p,2);
res = false(numPoints,1);
lb = infimum(I);
ub = supremum(I);

for iPoint = 1:numPoints
    p_curr = p(:,iPoint);
    res(iPoint) = all(lb <= p_curr) && all(ub >= p_curr);
end

%------------- END OF CODE --------------
