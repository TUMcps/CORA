function res = test_conPolyZono_polygon
% test_conPolyZono_polygon - unit test function of polygon
%
% Syntax:
%    res = test_conPolyZono_polygon
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
cPZ = conPolyZono.empty(2);
pgon = polygon(cPZ);
assert(representsa(pgon,"emptySet"));

% init cases
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];

% center
cPZ = conPolyZono(c);
pgon = polygon(cPZ);
[res,p] = representsa(pgon,'point');
assert(res & all(withinTol(p,c)));

% full set
cPZ = conPolyZono(c,G,E,A,b,EC);
pgon = polygon(cPZ);

% conversion is outer-approximative
xs = [cPZ.randPoint(100) cPZ.randPoint(100,'extreme')];
assert(all(pgon.contains(xs,'exact',1e-8)));

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
