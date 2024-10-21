function res = test_nonlinearReset_synchronize
% test_nonlinearReset_synchronize - test function for the synchronization
%    of a list of nonlinear reset functions
%
% Syntax:
%    res = test_nonlinearReset_synchronize
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
% See also: test_linearReset_synchronize

% Authors:       Mark Wetzlinger
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init nonlinear reset functions
f1 = @(x,u) [x(1)^2; -x(2)*u(1)];
nonlinReset1 = nonlinearReset(f1);
f2 = @(x,u) [x(1)*x(2) + u(1); -x(2)^2*u(1)];
nonlinReset2 = nonlinearReset(f2);

% lift for synchronize
n_high = 5; m_high = 3;
nonlinReset1 = lift(nonlinReset1,n_high,m_high,[1,2],1,false);
nonlinReset2 = lift(nonlinReset2,n_high,m_high,[4,5],3,false);

% synchronize
nonlinResets = [nonlinReset1;nonlinReset2];
nonlinReset_sync = synchronize(nonlinResets,[false;false;true;false;false]);
f_true = @(x,u) [x(1)^2; -x(2)*u(1); x(3); x(4)*x(5) + u(3); -x(5)^2*u(3)];
assert(isequalFunctionHandle(nonlinReset_sync.f,f_true));
assert(nonlinReset_sync.preStateDim == n_high);
assert(nonlinReset_sync.inputDim == m_high);
assert(nonlinReset_sync.postStateDim == n_high);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
