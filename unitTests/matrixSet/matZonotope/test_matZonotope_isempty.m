function res = test_matZonotope_isempty
% test_matZonotope_isempty - unit test function for emptiness checks
% 
% Syntax:
%    res = test_matZonotope_isempty
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
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

assert(isempty(matZonotope()));

% instantiate matrix zonotope
C = [0 0; 0 0];
G(:,:,1) = [1 3;1 2];
G(:,:,2) = [2 0; 1 -1];
matZ = matZonotope(C,G);

assert(~isempty(matZ));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
