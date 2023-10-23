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

resvec = [];

% empty
printMatrix([])
resvec(end+1) = true;

% scalar
printMatrix(1)
printMatrix(2)
printMatrix(-4)
printMatrix(Inf)
resvec(end+1) = true;

% row vector
printMatrix([1 2 3])
printMatrix([1 2 3 4 5])
resvec(end+1) = true;

% column vector
printMatrix([1 2 3]')
printMatrix([1 2 3 4 5]')
resvec(end+1) = true;

% matrix
printMatrix([1 2 3; 4 5 6])
resvec(end+1) = true;

% parameters
M = [1 2 3; 4 5 6];
printMatrix(M);
printMatrix(M,'%4.3f%s');
printMatrix(M,'high');
printMatrix(M,'high',true);
printMatrix(M,'high',false);
fprintf('\n')
resvec(end+1) = true;

% example
M = [2 3; -2 1];
printMatrix(M)
resvec(end+1) = true;

% combine results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
