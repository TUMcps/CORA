function res = test_nonlinearReset_lift
% test_nonlinearReset_lift - test function for the projection of a
%    nonlinear reset function
%
% Syntax:
%    res = test_nonlinearReset_lift
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
% Written:       13-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% only states
f = @(x,u) [x(1)*x(2); x(2)];
nonlinReset = nonlinearReset(f);
n_high = 6; m_high = 5; stateBind = [2,3]; inputBind = 3;
% ...with identity mapping
nonlinReset_lift = lift(nonlinReset,n_high,m_high,stateBind,inputBind,true);
f_true = @(x,u) [x(1); x(2)*x(3); x(3); x(4); x(5); x(6)];
assert(nonlinReset_lift.preStateDim == n_high);
assert(nonlinReset_lift.inputDim == m_high);
assert(nonlinReset_lift.postStateDim == n_high);
assert(isequalFunctionHandle(nonlinReset_lift.f,f_true));
% ...without identity mapping
nonlinReset_lift = lift(nonlinReset,n_high,m_high,stateBind,inputBind,false);
f_true = @(x,u) [0; x(2)*x(3); x(3); 0; 0; 0];
assert(nonlinReset_lift.preStateDim == n_high);
assert(nonlinReset_lift.inputDim == m_high);
assert(nonlinReset_lift.postStateDim == n_high);
assert(isequalFunctionHandle(nonlinReset_lift.f,f_true));

% states and inputs
f = @(x,u) [x(1) - u(1); x(1)*x(2)];
nonlinReset = nonlinearReset(f);
n_high = 6; m_high = 5; stateBind = [2,3]; inputBind = 3; id = true;
nonlinReset_lift = lift(nonlinReset,n_high,m_high,stateBind,inputBind,id);
f_true = @(x,u) [x(1); x(2) - u(3); x(2)*x(3); x(4); x(5); x(6)];
assert(nonlinReset_lift.preStateDim == n_high);
assert(nonlinReset_lift.inputDim == m_high);
assert(nonlinReset_lift.postStateDim == n_high);
assert(isequalFunctionHandle(nonlinReset_lift.f,f_true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
