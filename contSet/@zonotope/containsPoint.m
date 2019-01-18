function isContained = containsPoint(zonotope, p)
% containsPoint - determines if the point p is inside the zonotope
%
% As an optimization results is used, boundary points might not be
% recognized correctly.
%
% Syntax:  
%    result = containsPoint(zonotope, p)
%
% Inputs:
%    Z1 - zonotope object
%    p - point specified as a vector
%
% Outputs:
%    result - 1/0 if point is inside the zonotope or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

robustness = containsPointWithRobustness(zonotope, p, 0);
isContained = robustness >= 0;

