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

res = true;

% check empty conZonotope
cZ = conZonotope();
if dim(cZ) ~= 0
    res = false;
end

% simple case
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);

if dim(cZ) ~= 2
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
