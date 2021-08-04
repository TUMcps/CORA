function E = cartProd(E1,E2)
% dim - returns an overapproximation fo the cartesian product 
% between two ellipsoids
%
% Syntax:  
%    d = cartProd(E1,E2)
%
% Inputs:
%    E1,E2 - ellipsoid object
%
% Outputs:
%    E - ellipsoid object (result of cartesian product)
%
% Example: 
%    E1 = ellipsoid.generateRandom(2);
%    E2 = ellipsoid.generateRandom(2);
%    E = cartProd(E1,E2); 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Victor Gassmann
% Written:      19-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
if ~isa(E1,'ellipsoid') || ~isa(E2,'ellipsoid')
    error('Input arguments both need to be of type "ellipsoid"!');
end
IntE = cartProd(interval(E1),interval(E2));
r = rad(IntE);
q = center(IntE);
TOL = min(E1.TOL,E2.TOL);
% extract degenerate dimensions
ind_d = withinTol(r,zeros(size(r)),TOL);
n_nd = sum(~ind_d);

% construct final ellipsoid
E = ellipsoid(n_nd*diag(r.^2),q);
%------------- END OF CODE --------------