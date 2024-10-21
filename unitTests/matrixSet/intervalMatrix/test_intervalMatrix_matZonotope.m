function res = test_intervalMatrix_matZonotope
% test_intervalMatrix_matZonotope - unit test function for conversion to
%    matrix zonotopes 
% 
% Syntax:
%    res = test_intervalMatrix_matZonotope
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
% Written:       18-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% only one column
c = [1 2 3 4]';
d = [1 1 1 1]';
intMat = intervalMatrix(c,d);
matZ = matZonotope(intMat);
assert(all(dim(matZ) == [4,1]));
assert(matZ.numgens() == 4)
assert(compareMatrices(c,matZ.C))
assert(compareMatrices(matZ.G(:,:,1),[1;0;0;0]))
assert(compareMatrices(matZ.G(:,:,2),[0;1;0;0]))
assert(compareMatrices(matZ.G(:,:,3),[0;0;1;0]))
assert(compareMatrices(matZ.G(:,:,4),[0;0;0;1]));

% multiple rows and columns
c = [2 3 4; 5 6 7];
d = [1 0 1; 0 0 1];
intMat = intervalMatrix(c,d);
matZ = matZonotope(intMat);
assert(all(dim(matZ) == [2,3]));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
