function res = test_nonlinearReset_isequal
% test_nonlinearReset_isequal - test function for equality check of
%    nonlinearReset objects
%
% Syntax:
%    res = test_nonlinearReset_isequal
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
% Written:       09-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty case
nonlinReset_empty = nonlinearReset();
assert(isequal(nonlinReset_empty,nonlinReset_empty));

% only states
f = @(x,u) [x(1)*x(2); x(2)];
g = @(x,u) [x(2)*x(1); 1*x(2)];
nonlinReset = nonlinearReset(f);
nonlinReset_ = nonlinearReset(g);
assert(isequal(nonlinReset,nonlinReset_));

% % states and inputs
% f = @(x,u) [x(1) - u(1); x(1)*x(2)];
% nonlinReset = nonlinearReset(f);
% assert(nonlinReset.preStateDim == 2);
% assert(nonlinReset.inputDim == 1);
% assert(nonlinReset.postStateDim == 2);
% 
% % states and inputs, different output dimension
% f = @(x,u) x(1)*x(2) - u(1);
% nonlinReset = nonlinearReset(f);
% assert(nonlinReset.preStateDim == 2);
% assert(nonlinReset.inputDim == 1);
% assert(nonlinReset.postStateDim == 1);


% comparison nonlinearReset to linearReset
f_lin = @(x,u) [x(1) - u(1); x(2) + 2];
f_nonlin = @(x,u) [x(1) - u(1); x(1)*x(2) + 2];
nonlinReset_lin = nonlinearReset(f_lin);
nonlinReset_nonlin = nonlinearReset(f_nonlin);
A = [1 0; 0 1]; B = [-1; 0]; c = [0; 2];
linReset = linearReset(A,B,c);
assert(isequal(nonlinReset_lin,linReset));
assert(isequal(linReset,nonlinReset_lin));
assert(~isequal(nonlinReset_nonlin,linReset));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
