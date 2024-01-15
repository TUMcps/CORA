function res = test_intervalMatrix_contains
% test_intervalMatrix_contains - unit test function for contains
% 
% 
% Syntax:
%    res = test_intervalMatrix_contains
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

% instantiate interval matrix
c = [2 3 4; 5 6 0];
d = [1 0 1; 0 0 1];
intMat = intervalMatrix(c,d);

% center has to be contained
res(end+1,1) = contains(intMat,c);

% sample matrices
M = [0 3 5; 5 6 2];
res(end+1,1) = ~contains(intMat,M);
M = [2.5 3 3.5; 5 6 0.5];
res(end+1,1) = contains(intMat,M);

% include tolerance
M = [1 3 5; 5 6 1.1];
res(end+1,1) = ~contains(intMat,M,'exact',0.05);
res(end+1,1) = contains(intMat,M,'exact',0.2);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
