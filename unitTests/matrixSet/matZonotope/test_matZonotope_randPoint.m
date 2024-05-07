function res = test_matZonotope_randPoint
% test_matZonotope_randPoint - unit test function for random sampling
% 
% Syntax:
%    res = test_matZonotope_randPoint
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

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       03-April-2023
% Last update:   25-April-2024 (TL, better checks)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty matrix zonotope
matZ = matZonotope();
M = randPoint(matZ);
res = isnumeric(M) && isempty(M);

% instantiate matrix zonotope
C = [0 2; 1 -1; 1 -2];
G = []; G(:,:,1) = [1 1; -1 0; -2 1]; G(:,:,2) = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);

% try different syntaxes for runtime issues
M = randPoint(matZ);
res(end+1) = size(M,3) == 1;
M = randPoint(matZ,5);
res(end+1) = size(M,3) == 5;
M = randPoint(matZ,5,'extreme');
res(end+1) = size(M,3) == 5;
M = randPoint(matZ,5,'standard');
res(end+1) = size(M,3) == 5;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
