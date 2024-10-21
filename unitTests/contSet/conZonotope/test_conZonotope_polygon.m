function res = test_conZonotope_polygon
% test_conZonotope_polygon - unit test function of polygon
%
% Syntax:
%    res = test_conZonotope_polygon
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
cZ = conZonotope.empty(2);
pgon = polygon(cZ);
assert(representsa(pgon,"emptySet"));

% init cases
c = [0;1];
G = [3 0 1; 0 2 1];
A = [1 0 1]; b = 1;

% center
cZ = conZonotope(c);
pgon = polygon(cZ);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
cZ = conZonotope(c,G,A,b);
pgon = polygon(cZ);

% conversion is outer-approximative
xs = [cZ.randPoint(100) cZ.randPoint(100,'extreme')];
assert(all(pgon.contains(xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
