function res = test_ellipsoid_polygon
% test_ellipsoid_polygon - unit test function of polygon
%
% Syntax:
%    res = test_ellipsoid_polygon
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
E = ellipsoid.empty(2);
pgon = polygon(E);
assert(representsa(pgon,"emptySet"));

% init cases
Q = [ 5.4387811500952807 12.4977183618314545 ; 12.4977183618314545 29.6662117284481646 ];
q = [ -0.7445068341257537 ; 3.5800647524843665 ];

% center
E = ellipsoid(Q*0,q);
pgon = polygon(E);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,q,E.TOL)));

% full set
E = ellipsoid(Q, q, 0.000001);
pgon = polygon(E);

% conversion is inner-approximative
assert(E.contains(polytope(pgon)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
