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

res = true;

% instantiate interval matrix
intMat = intervalMatrix([2 3; 1 2],[1 0; 1 1]);
assert(~isempty(intMat));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
