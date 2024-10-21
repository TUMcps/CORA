function res = test_linearReset_lift
% test_linearReset_lift - test function for the projection of a linear
%    reset function
%
% Syntax:
%    res = test_linearReset_lift
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

% lift
n_high = 6; stateBind = [2,3];
m_high = 5; m_high_noB = 3;
inputBind = [2,3,4]; inputBind_noB = 2;
id = true;
linReset_A_lift = lift(linReset_A,n_high,m_high_noB,stateBind,inputBind_noB,id);
linReset_AB_lift = lift(linReset_AB,n_high,m_high,stateBind,inputBind,id);
linReset_ABc_lift = lift(linReset_ABc,n_high,m_high,stateBind,inputBind,id);

% check projected dimensions
assert(all(withinTol(linReset_A_lift.A(stateBind,stateBind),A,tol),'all'));
assert(all(withinTol(linReset_AB_lift.A(stateBind,stateBind),A,tol),'all'));
assert(all(withinTol(linReset_ABc_lift.A(stateBind,stateBind),A,tol),'all'));

% check new dimensions
otherDims = setdiff(1:n_high,stateBind); n_plus = n_high - numel(stateBind);
assert(all(withinTol(linReset_A_lift.A(otherDims,otherDims),eye(n_plus),tol),'all'));
assert(all(withinTol(linReset_AB_lift.A(otherDims,otherDims),eye(n_plus),tol),'all'));
assert(all(withinTol(linReset_ABc_lift.A(otherDims,otherDims),eye(n_plus),tol),'all'));

% check input matrix
assert(all(withinTol(linReset_A_lift.B,zeros(n_high,1),tol),'all'));
assert(all(withinTol(linReset_AB_lift.B(stateBind,inputBind),B,tol),'all'));
assert(all(withinTol(linReset_ABc_lift.B(stateBind,inputBind),B,tol),'all'));

% check vector
assert(all(withinTol(linReset_A_lift.c,zeros(n_high,1),tol)));
assert(all(withinTol(linReset_AB_lift.c,zeros(n_high,1),tol)));
assert(all(withinTol(linReset_ABc_lift.c(stateBind),c,tol)));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
