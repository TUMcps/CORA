function res = test_polytope_not
% test_polytope_not - unit test function of complement operation
%
% Syntax:
%    res = test_polytope_not
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
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true(0);

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
P_not = ~P;
res(end+1,1) = representsa(P_not,'emptySet');

% 1D, empty
A = [1; -1]; b = [1; -2];
P = polytope(A,b);
P_not = ~P;
res(end+1,1) = representsa(P_not,'fullspace');

% 2D, unbounded
A = [1 0]; b = 1;
P = polytope(A,b);
P_not = ~P;
P_true = polytope([-1 0],-1);
res(end+1,1) = isequal(P_not,P_true);

% 3D, fullspace
A = [0 0 0]; b = 2;
P = polytope(A,b);
P_not = ~P;
res(end+1,1) = representsa(P_not,'emptySet');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
