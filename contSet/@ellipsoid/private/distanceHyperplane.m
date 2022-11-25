function val = distanceHyperplane(E,H)
% distanceHyperplane - computes the distance from an ellipsoid to a
% conHyperplane object
%
% Syntax:  
%    res = distanceConHyperplane(E,H)
%
% Inputs:
%    E - ellipsoid object
%    H - conHyperplane object
%
% Outputs:
%    val - distance between E and H
%
% Example: 
%    ---
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal toolbox (ET).
% In Proceedings of the 45th IEEE Conference on Decision and Control (pp. 1498-1503). IEEE.
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
% can be <0
y = H.h.d; c = H.h.c;
q = E.q; Q = E.Q;
val = (abs(y-c'*q)-sqrt(c'*Q*c)) / sqrt(c'*c);
%------------- END OF CODE --------------