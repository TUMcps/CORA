function res = test_matZonotope_zonotope
% test_matZonotope_zonotope - unit test function for transpose
% 
% Syntax:
%    res = test_matZonotope_zonotope
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

% Authors:       Tobias Ladner
% Written:       20-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init matZonotope
C = [0 0; 0 0];
G(:,:,1) = [1 3; -1 2]; G(:,:,2) = [2 0; 1 -1];
matZ = matZonotope(C,G);

% check if center is contained
assert(contains(matZ,C));

% check if generators are contained (as center is origin)
assert(all(contains(matZ,G)));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
