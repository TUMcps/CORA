function res = test_capsule_ne
% test_capsule_ne - unit test function of '~='
%
% Syntax:
%    res = test_capsule_ne
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

% Authors:       Mark Wetzlinger
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);
tol = 1e-9;

% empty capsule
C = capsule.empty(2);
res(end+1,1) = ~(C ~= C);

% define properties
c = [2; 0; -1];
g = [0.2; -0.7; 0.4];
r = 1;
C = capsule(c,g,r);
C_tol = capsule(c,g,r+tol/2);
C_red = capsule(c(1:end-1),g(1:end-1),r);

% compare
res(end+1,1) = ~(C ~= C);
res(end+1,1) = ~ne(C,C);
res(end+1,1) = C ~= C_tol;
res(end+1,1) = ~ne(C,C_tol,tol);
res(end+1,1) = C ~= C_red;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
