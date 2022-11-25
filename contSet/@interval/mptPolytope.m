function poly = mptPolytope(obj)
% mptPolytope - convert an interval to a mptPolytope
%
% Syntax:  
%    poly = mptPolytope(obj)
%
% Inputs:
%    obj - interval object
%
% Outputs:
%    poly - mptPolytope object
%
% Example: 
%    int = interval([-1;0],[1;3]);
%    poly = mptPolytope(int);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/mptPolytope

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% check if the interval has the correct format
if size(obj,2) > 1
   error('Input argument ''obj'' has the wrong format!'); 
end

% construct halfspace constraints
n = length(obj);
C = [eye(n);-eye(n)];
d = [supremum(obj);-infimum(obj)];

% construct mptPolytope object
poly = mptPolytope(C,d);

%------------- END OF CODE --------------