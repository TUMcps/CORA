function res = test_conZonotope_copy
% test_conZonotope_copy - unit test function of copy
%
% Syntax:
%    res = test_conZonotope_copy
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
% Written:       02-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D conZonotope
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);
cZ_copy = copy(cZ);
assert(isequal(cZ,cZ_copy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
