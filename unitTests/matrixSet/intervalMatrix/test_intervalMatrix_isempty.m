function res = test_intervalMatrix_isempty
% test_intervalMatrix_isempty - unit test function for emptiness checks
% 
% Syntax:
%    res = test_intervalMatrix_isempty
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

% instantiate interval matrix
intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
res(end+1,1) = ~isempty(intMat);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
