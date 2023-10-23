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

% Authors:       Mark Wetzlinger
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty matrix zonotope
matZ = matZonotope();
M = randPoint(matZ);
res = iscell(M) && length(M) == 1 && isempty(M{1});

% instantiate matrix zonotope
C = [0 2; 1 -1; 1 -2];
G{1} = [1 1; -1 0; -2 1]; G{2} = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);

% try different syntaxes for runtime issues
randPoint(matZ);
randPoint(matZ,5);
randPoint(matZ,5,'extreme');
randPoint(matZ,5,'standard');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
