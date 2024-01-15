function res = test_conZonotope_dim
% test_conZonotope_dim - unit test function of dim
%
% Syntax:
%    res = test_conZonotope_dim
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
% Written:       14-March-2021
% Last update:   03-January-2023 (MW, moved random tests to testLongDuration)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
cZ = conZonotope.empty(3);
res(end+1,1) = dim(cZ) == 3;
cZ = conZonotope(zeros(2,0));
res(end+1,1) = dim(cZ) == 2;

% simple case
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
res(end+1,1) = dim(cZ) == 2;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
