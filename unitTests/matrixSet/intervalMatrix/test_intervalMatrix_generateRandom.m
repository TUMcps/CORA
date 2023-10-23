function res = test_intervalMatrix_generateRandom
% test_intervalMatrix_generateRandom - unit test function for random
%    generation of an interval matrix
% 
% 
% Syntax:
%    res = test_intervalMatrix_generateRandom
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

% no input arguments
intMat = intervalMatrix.generateRandom();

% dimension given
dims = [2,3];
intMat = intervalMatrix.generateRandom('Dimension',dims);
res = all(dim(intMat) == dims);

% center given
c = [2 3 4; 5 6 7];
intMat = intervalMatrix.generateRandom('Center',c);
res(end+1,1) = all(all(withinTol(center(intMat.int),c)));

% max radius given
r = 1;
intMat = intervalMatrix.generateRandom('MaxRadius',r);
res(end+1,1) = all(all(rad(intMat.int) < r)) | all(all(withinTol(rad(intMat.int),r,1e-10)));

% dimension and center given
intMat = intervalMatrix.generateRandom('Dimension',dims,'Center',c);
res(end+1,1) = all(dim(intMat) == dims) && all(all(withinTol(center(intMat.int),c)));

% center and max radius given
intMat = intervalMatrix.generateRandom('Center',c,'MaxRadius',r);
res(end+1,1) = all(all(withinTol(center(intMat.int),c))) && ...
    ( all(all(rad(intMat.int) < r)) | all(all(withinTol(rad(intMat.int),r,1e-10))) );


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
