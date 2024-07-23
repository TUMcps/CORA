function res = test_zonotope_reduceUnderApprox
% test_zonotope_reduceUnderApprox - unit test function of reduction
%    operation returning an inner approximation of the original zonotope
%
% Syntax:
%    res = test_zonotope_reduceUnderApprox
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
% Written:       20-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
% note: well copy calls (see below) to easier find specific wrong results

% empty zonotope: check all methods
Z = zonotope(zeros(2,0));
res(end+1,1) = representsa(reduceUnderApprox(Z,'sum',1),'emptySet');
res(end+1,1) = representsa(reduceUnderApprox(Z,'scale',1),'emptySet');
res(end+1,1) = representsa(reduceUnderApprox(Z,'linProg',1),'emptySet');
res(end+1,1) = representsa(reduceUnderApprox(Z,'wetzlinger',1),'emptySet');

% 2D
Z = zonotope([1;0],[1,3,-2,3,-1,0,2,3; ...
                    2,-1,1,0,3,-2,1,2]);
Z_red = reduceUnderApprox(Z,'sum',1);
res(end+1,1) = size(Z_red.G,2) == 2;
res(end+1,1) = contains(Z,Z_red,'exact',1e-10);
Z_red = reduceUnderApprox(Z,'scale',1);
res(end+1,1) = size(Z_red.G,2) == 2;
res(end+1,1) = contains(Z,Z_red,'exact',1e-10);
Z_red = reduceUnderApprox(Z,'linProg',1);
res(end+1,1) = size(Z_red.G,2) == 2;
res(end+1,1) = contains(Z,Z_red,'exact',1e-10);
Z_red = reduceUnderApprox(Z,'wetzlinger',1);
res(end+1,1) = size(Z_red.G,2) == 2;
res(end+1,1) = contains(Z,Z_red,'exact',1e-10);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
