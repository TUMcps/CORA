function res = test_zonotope_norm
% test_zonotope_norm - unit test function of norm
%
% Syntax:
%    res = test_zonotope_norm
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

% Authors:       Mark Wetzlinger, Victor Gassmann
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

TOL = 1e-6;

% empty case
Z_empty = zonotope.empty(2);
res = norm(Z_empty) == -Inf;

% full-dimensional case
c = zeros(2,1);
G = [2 5 4 3; -4 -6 2 3];
Z = zonotope(c,G);


% 2-norm test
val2_exact = norm(Z,2,'exact');
val2_ub = norm(Z,2,'ub');
val2_ubc = norm(Z,2,'ub_convex');

% compute vertices
V = vertices(Z);

% check exact vs. upper bound
if val2_exact > val2_ub && ~withinTol(val2_exact,val2_ub,TOL)
    res = false;
end

% check exact vs. upper bound (convex)
if val2_exact > val2_ubc && ~withinTol(val2_exact,val2_ubc,TOL)
    res = false;
end

% check exact vs. norm of all vertices
temp = abs(val2_exact-max(sqrt(sum(V.^2))))/val2_exact;
if temp > 0 && ~withinTol(temp,0,TOL)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
