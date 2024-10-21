function res = test_zonoBundle_polygon
% test_zonoBundle_polygon - unit test function of polygon
%
% Syntax:
%    res = test_zonoBundle_polygon
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
zB = zonoBundle.empty(2);
pgon = polygon(zB);
assert(representsa(pgon,"emptySet"));

% init cases
c = [1;2];

% center
zB = zonoBundle({zonotope(c)});
pgon = polygon(zB);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
Z1 = zonotope([1 3 0; 1 0 2]);
Z2 = zonotope([0 2 2; 0 2 -2]);
zB = zonoBundle({Z1,Z2});
pgon = polygon(zB);

% conversion is outer-approximative
xs = [zB.randPoint(100) zB.randPoint(100,'extreme')];
assert(all(pgon.contains(xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
