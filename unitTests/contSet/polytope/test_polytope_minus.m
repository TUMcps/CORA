function res = test_polytope_minus
% test_polytope_minus - unit test function of minus
%
% Syntax:
%    res = test_polytope_minus
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
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% 2D, vertex instantiation
V = [3 2; 0 3; -3 0; -1 -2; 2 -2]';
P = polytope(V);
% translate by vector
z = [2; 1];
P_minus = P - z;
V_true = V - z;
P_true = polytope(V_true);
res(end+1,1) = P_true == P_minus;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
