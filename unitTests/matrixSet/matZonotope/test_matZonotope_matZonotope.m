function res = test_matZonotope_matZonotope
% test_matZonotope_matZonotope - unit test function for constructor
% 
% Syntax:
%    res = test_matZonotope_matZonotope
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
% Written:       03-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% only center
C = 0;
matZ = matZonotope(C);
res(end+1,1) = compareMatrices(matZ.C,C);
C = [1; 1; 0];
matZ = matZonotope(C);
res(end+1,1) = compareMatrices(matZ.C,C);

% center and one generator
G = [];
G(:,:,1) = [1; 0; 0];
matZ = matZonotope(C,G);
res(end+1,1) = compareMatrices(matZ.C,C) && compareMatrices(matZ.G,G);

% center and multiple generators
C = [1 2 1; 3 2 0];
G = [];
G(:,:,1) = [2 0 1; -1 1 -2];
G(:,:,2) = [3 1 0; -1 -1 4];
G(:,:,3) = [0 1 -1; 3 1 2];
matZ = matZonotope(C,G);
res(end+1,1) = compareMatrices(matZ.C,C) && compareMatrices(matZ.G,G);

% copy constructor
matZ_ = matZonotope(matZ);
res(end+1,1) = true;

% conversion from zonotope: empty, only center, one/multiple generator(s)
Z = zonotope.empty(2);
matZ = matZonotope(Z);
Z = zonotope([1;2;1]);
matZ = matZonotope(Z);
Z = zonotope([1; 2],[1; 0]);
matZ = matZonotope(Z);
Z = zonotope([1; 2],[1 2 4 2; 0 2 -3 1]);
matZ = matZonotope(Z);

% legacy --- 

C = [1 2 1; 3 2 0];
G = [];
G(:,:,1) = [2 0 1; -1 1 -2];
G(:,:,2) = [3 1 0; -1 -1 4];
G(:,:,3) = [0 1 -1; 3 1 2];

G_legacy{1} = G(:,:,1);
G_legacy{2} = G(:,:,2);
G_legacy{3} = G(:,:,3);
matZ = matZonotope(C,G_legacy);
res(end+1,1) = compareMatrices(matZ.C,C) && compareMatrices(matZ.G,G);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
