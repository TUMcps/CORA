function res = test_printMatrix
% test_printMatrix - unit test function for printMatrix
%
% Syntax:
%    res = test_printMatrix()
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
% Written:       31-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty
printMatrix([])

% scalar
printMatrix(1)
printMatrix(2)
printMatrix(-4)
printMatrix(Inf)

% row vector
printMatrix([1 2 3])
printMatrix([1 2 3 4 5])

% column vector
printMatrix([1 2 3]')
printMatrix([1 2 3 4 5]')

% matrix
printMatrix([1 2 3; 4 5 6])

% parameters
M = [1 2 3; 4 5 6];
printMatrix(M);
printMatrix(M,'%.4e');
printMatrix(M,'high');
printMatrix(M,'high',true);
printMatrix(M,'high',false);
printMatrix(M,'high',true,false);
fprintf('\n')
printMatrix(M,'high',false,false);
fprintf('\n')

% example
M = [2 3; -2 1];
printMatrix(M)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
