function res = test_zonotope_polygon
% test_zonotope_polygon - unit test function of polygon
%
% Syntax:
%    res = test_zonotope_polygon
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
Z = zonotope.empty(2);
pgon = polygon(Z);
assert(representsa(pgon,"emptySet"));

% init cases
c = [1;2];

% center
Z = zonotope(c);
pgon = polygon(Z);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
Z = zonotope([1 3 0; 1 0 2]);
pgon = polygon(Z);

% conversion is outer-approximative
xs = [Z.randPoint(100) Z.randPoint(100,'extreme')];
assert(all(pgon.contains(xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
