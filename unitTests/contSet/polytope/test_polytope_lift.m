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

res = [];

% init polytope
A = [1 0;-1 0;0 1;0 -1;1 1];
b = [1;1;1;1;1];
P = polytope(A,b);

% lift
P_ = lift(P,10,[4,5]);

% true result
A_true = [zeros(5,3), A, zeros(5,5)];
P_true = polytope(A_true,b);

% compare results
res(end+1,1) = P_ == P_true;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
