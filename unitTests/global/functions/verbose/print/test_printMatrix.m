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
printMatrix([]);

% scalar
printMatrix(1);
printMatrix(2);
printMatrix(-4);
printMatrix(Inf);

% row vector
printMatrix([1 2 3]);
printMatrix([1 2 3 4 5]);

% column vector
printMatrix([1 2 3]');
printMatrix([1 2 3 4 5]');

% matrix
printMatrix([1 2 3; 4 5 6]);

% non-integer
printMatrix(rand(2,3));

% special matrices
printMatrix(zeros(2,3));
printMatrix(ones(2,3));
printMatrix(4*ones(2,3));
printMatrix(eye(2));

% parameters
M = [1.1 2 3.3; 4 5 6];
printMatrix(M);
printMatrix(M,'%.4e');
printMatrix(M,'high');
printMatrix(M,'high',true);
printMatrix(M,'high',false);
printMatrix(M,'high',true,false);
fprintf('\n')
printMatrix(M,'high',false,false);
fprintf('\n')

% test fid
filename = 'test.txt';
printMatrix(filename,M,'high');
M_copy = eval(fileread(filename));
assert(isequal(M,M_copy));
delete(filename)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
