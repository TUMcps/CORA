function res = test_intervalMatrix_display
% test_intervalMatrix_display - unit test function for display (only check
%    for runtime errors) 
% 
% Syntax:
%    res = test_intervalMatrix_display
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

res = true;

% scalar
intMat = intervalMatrix(1,0)

% vector
c = [1; 0; 1]; d = [1; 2; 2];
intMat = intervalMatrix(c,d)

% matrix
c = [2 3 4; 5 6 0];
d = [1 0 1; 0 0 1];
intMat = intervalMatrix(c,d)

% ------------------------------ END OF CODE ------------------------------
