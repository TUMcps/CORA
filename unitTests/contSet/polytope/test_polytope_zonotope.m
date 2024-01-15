function res = test_polytope_zonotope
% test_polytope_zonotope - unit test function of zonotope conversion
%
% Syntax:
%    res = test_polytope_zonotope
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
% Written:       08-January-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
P = polytope.empty(2);
Z = zonotope(P);
res(end+1,1) = representsa(Z,'emptySet') && dim(Z) == 2;


% 2D, bounded
A = [1 0; -1 1; -1 -1]; b = [1;1;1];
P = polytope(A,b);
Z = zonotope(P);
res(end+1,1) = contains(Z,P);


% combine results
res = all(res);


% unsupported
try
    % unbounded
    P = polytope.Inf(2);
    Z = zonotope(P);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
