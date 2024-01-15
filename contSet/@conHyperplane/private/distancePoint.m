function val = distancePoint(hyp,p)
% distancePoint - computes the distance from a constrained hyperplane to a
%    single point
%
% Syntax:
%    res = distancePoint(hyp,p)
%
% Inputs:
%    hyp - conHyperplane object
%    p - point
%
% Outputs:
%    val - distance
%
% Example:
%    hyp = conHyperplane(halfspace([1;1],0),[1 0;-1 0],[2;2])
%    p = [1;1];
%    distance(hyp,p) % calls distancePoint
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% if conHyperplane actually is a hyperplane => analytical solution
if isempty(hyp.C) || all(all(hyp.C))
    val = abs(hyp.b-hyp.a*p)/(hyp.a*hyp.a);
else
    n = length(p);
    % solve optimization problem: ||x-y||_2^2, s.t. y\in{x|C*x<=d}
    val = quadprog(2*eye(n),-2*p,hyp.C,hyp.d);
end

% ------------------------------ END OF CODE ------------------------------
