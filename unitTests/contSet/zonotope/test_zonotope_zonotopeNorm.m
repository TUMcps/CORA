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

res = true(0);

% empty zonotopes
Z = zonotope.empty(2);
res(end+1,1) = zonotopeNorm(Z,[1;-1]) == Inf;
res(end+1,1) = zonotopeNorm(Z,double.empty(2,0)) == 0;

% 2D, only center
Z = zonotope([0;0]);
res(end+1,1) = zonotopeNorm(Z,[0;0]) == 0;
res(end+1,1) = zonotopeNorm(Z,[1;0]) == Inf;

% 2D, center and generators
c = [0;0]; G = [1 -2 2 0; -1 1 0 1];
Z = zonotope(c,G);
res(end+1,1) = withinTol(zonotopeNorm(Z,[5;3]),2.2);
res(end+1,1) = withinTol(zonotopeNorm(Z,[-5;3]),1);
% shifted center does not influence the result
Z = Z + [1;-2];
res(end+1,1) = withinTol(zonotopeNorm(Z,[5;3]),2.2);
res(end+1,1) = withinTol(zonotopeNorm(Z,[-5;3]),1);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
