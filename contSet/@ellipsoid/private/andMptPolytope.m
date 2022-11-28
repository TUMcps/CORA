function E = andMptPolytope(E,P,mode)
% andMptPolytope - Computes an inner-approximation or outer-approximation of
%    the intersection between an ellipsoid and a polytope
%
% Syntax:  
%    E = andMptPolytope(E,P,mode)
%
% Inputs:
%    E - ellipsoid object
%    P - mptPolytope object
%    mode - type of approximation
%               'inner'
%               'outer'
%
% Outputs:
%    E - ellipsoid after intersection
%
% References:
%   [1] Kurzhanskiy, A.A. and Varaiya, P., 2006, December. Ellipsoidal
%       toolbox (ET). In Proceedings of the 45th IEEE Conference o
%       Decision and Control (pp. 1498-1503). IEEE.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Victor Gassmann
% Written:      07-June-2022
% Last update:  05-July-2022 (VG: removed input checks; now in parent function)
% Last revision:---

%------------- BEGIN CODE --------------

% halfspace rep of P (is computed if not there; expensive!)
A = P.P.A;
b = P.P.b;

% loop over each halfspace and compute using ellipsoidHalfspace
for i=1:size(A,1)
    E = andHalfspace(E,halfspace(A(i,:)',b(i)),mode);
end

%------------- END OF CODE --------------