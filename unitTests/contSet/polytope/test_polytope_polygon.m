function res = test_polytope_polygon
% test_polytope_polygon - unit test function of polygon
%
% Syntax:
%    res = test_polytope_polygon
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
P = polytope.empty(2);
pgon = polygon(P);
assert(representsa(pgon,"emptySet"));

% init cases
A = [1 0 -1 0 1; 0 1 0 -1 1]';
b = [3; 2; 3; 2; 1];
c = [1;2];

% center
P = polytope(c);
pgon = polygon(P);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
P = polytope(A,b);
pgon = polygon(P);

% conversion is outer-approximative
xs = [P.randPoint(100) P.randPoint(100,'extreme')];
assert(all(pgon.contains(xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
