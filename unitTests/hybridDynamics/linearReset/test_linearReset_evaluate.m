function res = test_linearReset_evaluate
% test_linearReset_evaluate - test function for the evaluation of a linear
%    reset function
%
% Syntax:
%    res = test_linearReset_evaluate
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
% See also: test_nonlinearReset_evaluate

% Authors:       Mark Wetzlinger
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = eps;

% init linear reset functions
A = [1 2; 0 -1];
B = [2 0 1; -1 0 0];
c = [1; -5];
linReset_A = linearReset(A);
linReset_AB = linearReset(A,B);
linReset_ABc = linearReset(A,B,c);

% vectors
x = [5;-2];
u = [1; 4; -3];

x_ = evaluate(linReset_A,x);
x_true = A*x;
assert(all(withinTol(x_,x_true,tol)));
x_ = evaluate(linReset_AB,x);
assert(all(withinTol(x_,x_true,tol)));

x_ = evaluate(linReset_AB,x);
x_true = A*x;
assert(all(withinTol(x_,x_true,tol)));
x_ = evaluate(linReset_AB,x,u);
x_true = A*x + B*u;
assert(all(withinTol(x_,x_true,tol)));

x_ = evaluate(linReset_ABc,x);
x_true = A*x + c;
assert(all(withinTol(x_,x_true,tol)));
x_ = evaluate(linReset_ABc,x,u);
x_true = A*x + B*u + c;
assert(all(withinTol(x_,x_true,tol)));

% zonotopes
Z_x = zonotope([1;-1],[1 0 -1; 2 1 1]);
Z_u = zonotope([0;1;-1],[1 0 1; -1 -2 1; 2 1 3]);

x_ = evaluate(linReset_A,Z_x);
x_true = A*Z_x;
assert(isequal(x_,x_true,tol));
x_ = evaluate(linReset_AB,Z_x);
assert(isequal(x_,x_true,tol));

x_ = evaluate(linReset_AB,Z_x);
x_true = A*Z_x;
assert(isequal(x_,x_true,tol));
x_ = evaluate(linReset_AB,Z_x,Z_u);
x_true = A*Z_x + B*Z_u;
assert(isequal(x_,x_true,tol));

x_ = evaluate(linReset_ABc,Z_x);
x_true = A*Z_x + c;
assert(isequal(x_,x_true,tol));
x_ = evaluate(linReset_ABc,Z_x,Z_u);
x_true = A*Z_x + B*Z_u + c;
assert(isequal(x_,x_true,tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
