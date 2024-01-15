function res = test_polytope_lift
% test_polytope_lift - unit test function of lift
%
% Syntax:
%    res = test_polytope_lift
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
% Written:       03-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 1D, fully empty
A = zeros(0,1); b = zeros(0,0);
P = polytope(A,b);
P_lift = lift(P,3,2);
% true result
A_true = zeros(0,3); b_true = zeros(0,0);
P_true = polytope(A_true,b_true);
res(end+1,1) = P_lift == P_true;


% 2D, bounded
A = [1 0;-1 0;0 1;0 -1;1 1]; b = [1;1;1;1;1];
P = polytope(A,b);
P_lift = lift(P,10,[4,5]);
% true result
A_true = [zeros(5,3), A, zeros(5,5)]; b_true = b;
P_true = polytope(A_true,b_true);
res(end+1,1) = P_lift == P_true;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
