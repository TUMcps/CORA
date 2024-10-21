function res = test_intervalMatrix_printSet
% test_intervalMatrix_printSet - unit test function of printSet
%
% Syntax:
%    res = test_intervalMatrix_printSet
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

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test normal set
C = [0 2; 3 1];
D = [1 2; 1 1];
intMat = intervalMatrix(C,D);

printSet(intMat)
printSet(intMat,'high')
printSet(intMat,'high',true)
printSet(intMat,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
