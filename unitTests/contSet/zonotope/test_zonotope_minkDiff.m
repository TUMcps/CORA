function res = test_zonotope_minkDiff
% test_zonotope_minkDiff - unit test function of Minkowski difference
%
% Syntax:
%    res = test_zonotope_minkDiff
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
% Written:       05-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

tol = 1e-8;

% 2D, only centers
Z1 = zonotope([1;-1]);
Z2 = zonotope([0;1]);
Z_diff = minkDiff(Z1,Z2);
Z_true = zonotope([1;-2]);
assert(isequal(Z_diff,Z_true));

% 2D, box - 1/2*box
Z1 = zonotope(zeros(2,1),eye(2));
Z2 = zonotope(zeros(2,1),0.7*eye(2));
Z_diff = minkDiff(Z1,Z2);
Z_true = zonotope(zeros(2,1),0.3*eye(2));
assert(isequal(Z_diff,Z_true));

% 2D, aligned generators
G1 = [1 0; -1 2; 3 -1; 0 2]';
G2 = 0.7*G1;
G_diff = G1 - G2;
G2 = G2(:,[3, 2, 4, 1]);
Z1 = zonotope(zeros(2,1),G1);
Z2 = zonotope(zeros(2,1),G2);
Z_diff = minkDiff(Z1,Z2);
Z_true = zonotope(zeros(2,1),G_diff);
assert(isequal(Z_diff,Z_true,tol));

% 2D, different methods (Mink. diff of 2D zonotopes is closed!)
G1 = [1 0; -1 2; 3 -1; 0 2]';
Z1 = zonotope(zeros(2,1),G1);
G2 = 0.1*[1 1; -1 3; 5 -2; 1 2]';
Z2 = zonotope(zeros(2,1),G2);
P_true = minkDiff(polytope(Z1),polytope(Z2));
% ...inner approximations
Z_diff = minkDiff(Z1,Z2,'inner');
assert(isequal(P_true,Z_diff,tol));
Z_diff = minkDiff(Z1,Z2,'inner:conZonotope');
assert(isequal(P_true,Z_diff,tol));
Z_diff = minkDiff(Z1,Z2,'inner:RaghuramanKoeln');
assert(isequal(P_true,Z_diff,tol));
% ...outer approximations
Z_diff = minkDiff(Z1,Z2,'outer');
assert(isequal(P_true,Z_diff,tol));
Z_diff = minkDiff(Z1,Z2,'outer:coarse');
assert(isequal(P_true,Z_diff,tol));
Z_diff = minkDiff(Z1,Z2,'outer:scaling');
assert(isequal(P_true,Z_diff,tol));

% 2D, empty result
G1 = [1 0; -1 2; 3 -1; 0 2]';
Z1 = zonotope(zeros(2,1),G1);
G2 = [-1 0; 1 2; -3 -1; 0 2]';
Z2 = zonotope(zeros(2,1),G2);
Z_diff = minkDiff(Z1,Z2);
assert(representsa_(Z_diff,'emptySet',tol));
Z_diff = minkDiff(Z1,Z2,'inner');
assert(representsa_(Z_diff,'emptySet',tol));
Z_diff = minkDiff(Z1,Z2,'inner:conZonotope');
assert(representsa_(Z_diff,'emptySet',tol));
Z_diff = minkDiff(Z1,Z2,'inner:RaghuramanKoeln');
assert(representsa_(Z_diff,'emptySet',tol));
Z_diff = minkDiff(Z1,Z2,'outer');
assert(representsa_(Z_diff,'emptySet',tol));
Z_diff = minkDiff(Z1,Z2,'outer:coarse');
assert(representsa_(Z_diff,'emptySet',tol));
Z_diff = minkDiff(Z1,Z2,'outer:scaling');
assert(representsa_(Z_diff,'emptySet',tol));

% 2D, degenerate result
Z1 = zonotope(zeros(2,1),eye(2));
Z2 = zonotope(zeros(2,1),[1;1]);
Z_diff = minkDiff(Z1,Z2);
assert(representsa_(Z_diff,'origin',tol));
Z_diff = minkDiff(Z1,Z2,'inner');
assert(representsa_(Z_diff,'origin',tol));
Z_diff = minkDiff(Z1,Z2,'inner:conZonotope');
assert(representsa_(Z_diff,'origin',tol));
Z_diff = minkDiff(Z1,Z2,'inner:RaghuramanKoeln');
assert(representsa_(Z_diff,'origin',tol));
Z_diff = minkDiff(Z1,Z2,'outer');
assert(representsa_(Z_diff,'origin',tol));
Z_diff = minkDiff(Z1,Z2,'outer:coarse');
assert(representsa_(Z_diff,'origin',tol));
Z_diff = minkDiff(Z1,Z2,'outer:scaling');
assert(representsa_(Z_diff,'origin',tol));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
