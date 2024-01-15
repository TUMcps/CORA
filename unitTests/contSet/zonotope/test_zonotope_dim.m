function res = test_zonotope_dim
% test_zonotope_dim - unit test function of dim
%
% Syntax:
%    res = test_zonotope_dim
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
% Written:       27-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% check empty zonotope
Z = zonotope.empty(2);
res(end+1,1) = dim(Z) == 2;
Z = zonotope(zeros(3,0));
res(end+1,1) = dim(Z) == 3;


% instantiate zonotope
c = [-2; 1]; G = [2 4 5 3 3; 0 3 5 2 3];
Z = zonotope(c,G);
res(end+1,1) = dim(Z) == 2;

% no generator matrix
Z = zonotope(c);
res(end+1,1) = dim(Z) == 2;


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
