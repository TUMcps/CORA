function [result] = in(obj1,obj2)
% in - determines if elements of a zonotope obj2 intersect a constrained
%      hyperplane
%
% Syntax:  
%    [result] = in(obj1,obj2)
%
% Inputs:
%    obj1 - constrainedHyperplane object
%    obj2 - zonotope object
%
% Outputs:
%    result - 1/0 if zonotope is in, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      22-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parts of the zonotope are located in the constrained hyperplane if the
% zonotope intersects the hyperplane
result = isIntersecting(obj1,obj2);

%------------- END OF CODE --------------