function res = test_taylorLinSys_eAdt
% test_taylorLinSys_eAdt - unit test for helper class taylorLinSys
%
% Syntax:
%    res = test_taylorLinSys_eAdt
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       18-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% init matrix and object
A = [16 -2 3 13; 5 -11 10 -8; -9 7 6 12; -4 14 -15 1];
tls = taylorLinSys(A);

% init time step size
timeStep = 0.01;
options = struct('timeStep',timeStep);

% compute exponential matrix and compare to directly computed solution
val = computeField(tls,'eAdt',options);
val_true = expm(A*timeStep);
assert(compareMatrices(val,val_true,tol,"equal",true));

% different time step size
timeStep = 0.05;
options = struct('timeStep',timeStep);
val = computeField(tls,'eAdt',options);
val_true = expm(A*timeStep);
assert(compareMatrices(val,val_true,tol,"equal",true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
