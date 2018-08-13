function [result] = in(obj1,obj2)
% in - determines if elements of a zonotope obj2 are in an interval 
%
% Syntax:  
%    [result] = in(obj1,obj2)
%
% Inputs:
%    obj1 - interval object
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
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% obtain constraints C*x<d for the interval
n = length(obj1);
C = [eye(n);-eye(n)];
d = [supremum(obj1);-infimum(obj1)];

%affine map
Znew = C*obj2 + (-d);

%compute interval
IHnew = interval(Znew);

if all(infimum(IHnew)<=0)
    result = 1;
else
    result = 0;
end


%------------- END OF CODE --------------