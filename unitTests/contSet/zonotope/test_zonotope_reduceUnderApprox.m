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
% Last update:   06-March-2026 (LL, update reduction methods)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: we'll copy calls (see below) to easier find specific wrong results

% empty zonotope: check all methods
Z = zonotope(zeros(2,0));
assert(representsa(reduceUnderApprox(Z,'yang',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'raghuraman',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'kochdumper',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'sadraddini',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'cluster',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'boxCone',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'scale',1),'emptySet'));
assert(representsa(reduceUnderApprox(Z,'nlp',1),'emptySet'));

% 2D
Z = zonotope([1;0],[1,3,-2,3,-1,0,2,3; ...
                    2,-1,1,0,3,-2,1,2]);

% test method yang
Z_red = reduceUnderApprox(Z,'yang',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method raghuraman
Z_red = reduceUnderApprox(Z,'raghuraman',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method kochdumper
Z_red = reduceUnderApprox(Z,'kochdumper',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method sadraddini
Z_red = reduceUnderApprox(Z,'sadraddini',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method cluster
Z_red = reduceUnderApprox(Z,'cluster',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method boxCone
Z_red = reduceUnderApprox(Z,'boxCone',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method scale
Z_red = reduceUnderApprox(Z,'scale',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));

% test method nlp
Z_red = reduceUnderApprox(Z,'nlp',1);
assert(size(Z_red.G,2) == 2);
assert(contains(Z,Z_red,'exact',1e-10));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
