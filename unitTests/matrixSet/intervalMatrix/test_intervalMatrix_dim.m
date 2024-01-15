function res = test_intervalMatrix_dim
% test_intervalMatrix_dim - unit test function for dimension
% 
% 
% Syntax:
%    res = test_intervalMatrix_dim
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

res = true(0);

% only one row
c = [1 2 3 4];
d = [1 1 1 1];
intMat = intervalMatrix(c,d);
res(end+1,1) = all(dim(intMat) == [1,4]);

% only one column
c = [1 2 3 4]';
d = [1 1 1 1]';
intMat = intervalMatrix(c,d);
res(end+1,1) = all(dim(intMat) == [4,1]);


% multiple rows and columns
c = [2 3 4; 5 6 7];
d = [1 0 1; 0 0 1];
intMat = intervalMatrix(c,d);
res(end+1,1) = all(dim(intMat) == [2,3]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
