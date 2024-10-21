function res = test_zonotope_zonotopeNorm
% test_zonotope_zonotopeNorm - unit test function of zonotope norm
%
% Syntax:
%    res = test_zonotope_zonotopeNorm
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty zonotopes
Z = zonotope.empty(2);
assert(zonotopeNorm(Z,[1;-1]) == Inf);
assert(zonotopeNorm(Z,double.empty(2,0)) == 0);

% 2D, only center
Z = zonotope([0;0]);
assert(zonotopeNorm(Z,[0;0]) == 0);
assert(zonotopeNorm(Z,[1;0]) == Inf);

% 2D, center and generators
c = [0;0]; G = [1 -2 2 0; -1 1 0 1];
Z = zonotope(c,G);
assert(withinTol(zonotopeNorm(Z,[5;3]),2.2));
assert(withinTol(zonotopeNorm(Z,[-5;3]),1));
% shifted center does not influence the result
Z = Z + [1;-2];
assert(withinTol(zonotopeNorm(Z,[5;3]),2.2));
assert(withinTol(zonotopeNorm(Z,[-5;3]),1));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
