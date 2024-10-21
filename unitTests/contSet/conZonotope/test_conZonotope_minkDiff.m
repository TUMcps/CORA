function res = test_conZonotope_minkDiff
% test_conZonotope_minkDiff - unit test function for the computation of the
%    Minkowski difference between a constrained zonotope and another set
%
% Syntax:
%    res = test_conZonotope_minkDiff
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
% Written:       02-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-6;

% init constrained zonotope
c = [1; 0];
G = [1.5 -1.5 0.5;1 0.5 -1];
A = [1 1 1]; b = 1;
cZ = conZonotope(c,G,A,b);
P = polytope(cZ);

% Minkowski difference with vector
v = [2; -1];
cZ_minkDiff = minkDiff(cZ,v);
assert(all(withinTol(cZ.c - v, cZ_minkDiff.c,tol)));

% Minkowski difference with set
I = interval([-0.5;0],[1;2]);
cZ_minkDiff = minkDiff(cZ,I,'exact');
cZ_true = minkDiff(P,I);
assert(isequal(polytope(cZ_minkDiff),cZ_true,tol));

cZ_minkDiff = minkDiff(cZ,I,'inner:vertices');
cZ_true = minkDiff(P,I);
assert(contains(cZ_true,cZ_minkDiff,'exact',tol));

cZ_minkDiff = minkDiff(cZ,I,'inner:Vinod');
cZ_true = minkDiff(P,I);
assert(contains(cZ_true,cZ_minkDiff,'exact',tol));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
