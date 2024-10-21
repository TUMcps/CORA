function res = test_linearReset_isequal
% test_linearReset_isequal - test function for equality check of
%    linearReset objects
%
% Syntax:
%    res = test_linearReset_isequal
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
% See also: test_nonlinearReset_isequal

% Authors:       Mark Wetzlinger
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty cases
linReset_empty = linearReset();
assert(isequal(linReset_empty,linReset_empty));
linReset = linearReset([],[],[]);
assert(isequal(linReset,linReset_empty));
assert(isequal(linReset_empty,linReset));

% only state matrix
A = [1 2; 0 -1];
linReset = linearReset(A);
assert(isequal(linReset,linReset));
A_tol = [1+1e-10 2; 0 -1-1e-10];
linReset_tol = linearReset(A_tol);
assert(~isequal(linReset,linReset_tol));
assert(isequal(linReset,linReset_tol,1e-8));

% state and input matrix
A = [1 2; 0 -1];
B = [2 0 1; -1 0 0];
linReset = linearReset(A,B);
assert(isequal(linReset,linReset));
A_tol = [1+1e-10 2; 0 -1-1e-10];
B_tol = [2 0 1; -1 1e-10 0];
linReset_tol = linearReset(A_tol,B_tol);
assert(~isequal(linReset,linReset_tol));
assert(isequal(linReset,linReset_tol,1e-8));

% state matrix, input matrix, offset vector
A = [1 2; 0 -1];
B = [2 0 1; -1 0 0];
c = [-1; 0];
linReset = linearReset(A,B,c);
assert(isequal(linReset,linReset));
A_tol = [1+1e-10 2; 0 -1-1e-10];
B_tol = [2 0 1; -1 1e-10 0];
c_tol = [-1; 1e-10];
linReset_tol = linearReset(A_tol,B_tol,c_tol);
assert(~isequal(linReset,linReset_tol));
assert(isequal(linReset,linReset_tol,1e-8));

% default cases
A = [1 0; 0 -1];
B_def = [0; 0];
c_def = [0; 0];
linReset = linearReset(A);
linReset_def1 = linearReset(A,B_def);
linReset_def2 = linearReset(A,[],c_def);
linReset_def3 = linearReset(A,B_def,c_def);
assert(isequal(linReset,linReset_def1));
assert(isequal(linReset,linReset_def2));
assert(isequal(linReset,linReset_def3));

% all-zero
A = zeros(3); B = zeros(3,1); c = zeros(3,1);
linReset = linearReset(A,B,c);
A_tol = [0 0 1e-10; 0 0 0; -1e-10 0 0];
B_tol = [0; 0; 1e-10];
c_tol = [-1e-10; 0; 0];
linReset_tol = linearReset(A_tol,B_tol,c_tol);
assert(~isequal(linReset,linReset_tol));
assert(isequal(linReset,linReset_tol,1e-10));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
