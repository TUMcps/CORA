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

res = true(0);

% only one column
c = [1 2 3 4]';
d = [1 1 1 1]';
intMat = intervalMatrix(c,d);
matZ = matZonotope(intMat);
res(end+1,1) = all(dim(matZ) == [4,1]);
res(end+1,1) = matZ.gens == 4 && compareMatrices(c,matZ.center) ...
    && compareMatrices(matZ.generator{1},[1;0;0;0]) ...
    && compareMatrices(matZ.generator{2},[0;1;0;0]) ...
    && compareMatrices(matZ.generator{3},[0;0;1;0]) ...
    && compareMatrices(matZ.generator{4},[0;0;0;1]);

% multiple rows and columns
c = [2 3 4; 5 6 7];
d = [1 0 1; 0 0 1];
intMat = intervalMatrix(c,d);
matZ = matZonotope(intMat);
res(end+1,1) = all(dim(matZ) == [2,3]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
