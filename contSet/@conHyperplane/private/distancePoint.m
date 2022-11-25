function val = distancePoint(H,x)
% distancePoint - computes the distance from a conHyperplane to a point
%
% Syntax:  
%    res = distance(H,x)
%
% Inputs:
%    H - conHyperplane object
%    x - point
%
% Outputs:
%    val - distance between H and x
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      08-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
% if conHyperplane actually is a hyperplane => analytical solution
if isempty(H.C) || all(all(H.C))
    val = abs(H.h.d-H.h.c'*x)/(H.h.c'*H.h.c);
else
    n = length(x);
    % solve optimization problem: ||x-y||_2^2, s.t. y\in{x|C*x<=d}
    val = quadprog(2*eye(n),-2*x,H.C,H.d);
end
%------------- END OF CODE --------------