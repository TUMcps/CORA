function val = priv_distanceHyperplane(E,P)
% priv_distanceHyperplane - computes the distance from an ellipsoid to a
%    hyperplane
%
% Syntax:
%    res = priv_distanceHyperplane(E,P)
%
% Inputs:
%    E - ellipsoid object
%    P - polytope object representing a hyperplane
%
% Outputs:
%    val - distance between ellipsoid and hyperplane
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
%       In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
%   [2] https://www2.eecs.berkeley.edu/Pubs/TechRpts/2006/EECS-2006-46.pdf
%       (for details)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/distance

% Authors:       Victor Gassmann
% Written:       08-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% from [2, Sec. 2.1, Eq. (2.10)]
% can be <0
y = P.be; c = P.Ae';
q = E.q; Q = E.Q;
val = (abs(y-c'*q)-sqrt(c'*Q*c)) / sqrt(c'*c);

% ------------------------------ END OF CODE ------------------------------
