function res = test_matZonotope_dim
% test_matZonotope_dim - unit test function for dimension read
% 
% Syntax:
%    res = test_matZonotope_dim
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
res = all(dim(matZ) == [0,0]);

% scalar
C = 0;
G{1} = 1; G{2} = -2;
matZ = matZonotope(C,G);
res(end+1,1) = all(dim(matZ) == [1,1]);

% nx1 vector
C = [0; 1; 1];
G{1} = [1; -1; -2]; G{2} = [-2; 0; 1];
matZ = matZonotope(C,G);
res(end+1,1) = all(dim(matZ) == [3,1]);
res(end+1,1) = dim(matZ,1) == 3;

% matrix
C = [0 2; 1 -1; 1 -2];
G{1} = [1 1; -1 0; -2 1]; G{2} = [-2 0; 0 1; 1 -1];
matZ = matZonotope(C,G);
res(end+1,1) = all(dim(matZ) == [3,2]);
res(end+1,1) = dim(matZ,1) == 3;
res(end+1,1) = dim(matZ,2) == 2;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
