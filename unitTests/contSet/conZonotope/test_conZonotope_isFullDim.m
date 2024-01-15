function res = test_conZonotope_isFullDim
% test_conZonotope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = test_conZonotope_isFullDim
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
% Written:       21-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% check empty conZonotope object
cZ = conZonotope.empty(2);
res(end+1,1) = ~isFullDim(cZ);

% constrained zonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
res(end+1,1) = isFullDim(cZ);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
