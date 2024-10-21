function res = test_taylorLinSys_Apower_Apower_abs
% test_taylorLinSys_Apower_Apower_abs - unit test for helper class
%    taylorLinSys
%
% Syntax:
%    res = test_taylorLinSys_Apower_Apower_abs
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

% compute first power: A^1*dt^1/1!
timeStep = 0.05; ithpower = 1;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_dt_fact',options);
val_true = A * timeStep;
assert(compareMatrices(val,val_true,tol,"equal",true));

% compute next-higher power: A^2*dt^2/2!
timeStep = 0.05; ithpower = 2;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_dt_fact',options);
val_true = A^2 * timeStep^2 / 2;
assert(compareMatrices(val,val_true,tol,"equal",true));

% compute power with gap: A^5*dt^5/5!
timeStep = 0.05; ithpower = 5;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_dt_fact',options);
val_true = A^5 * timeStep^5 / factorial(5);
assert(compareMatrices(val,val_true,tol,"equal",true));
% read out same value again
val = computeField(tls,'Apower_dt_fact',options);
assert(compareMatrices(val,val_true,tol,"equal",true));

% different time size, start from ithpower > 1
timeStep = 0.1; ithpower = 3;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_dt_fact',options);
val_true = A^3 * timeStep^3 / factorial(3);
assert(compareMatrices(val,val_true,tol,"equal",true));


% same tests with absolute value of A
A_abs = abs(A);
% compute first power: |A|^1*dt^1/1!
timeStep = 0.05; ithpower = 1;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_abs_dt_fact',options);
val_true = A_abs * timeStep;
assert(compareMatrices(val,val_true,tol,"equal",true));

% compute next-higher power: |A|^2*dt^2/2!
timeStep = 0.05; ithpower = 2;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_abs_dt_fact',options);
val_true = A_abs^2 * timeStep^2 / 2;
assert(compareMatrices(val,val_true,tol,"equal",true));

% compute power with gap: |A|^5*dt^5/5!
timeStep = 0.05; ithpower = 5;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_abs_dt_fact',options);
val_true = A_abs^5 * timeStep^5 / factorial(5);
assert(compareMatrices(val,val_true,tol,"equal",true));
% read out same value again
val = computeField(tls,'Apower_abs_dt_fact',options);
assert(compareMatrices(val,val_true,tol,"equal",true));

% different time size, start from ithpower > 1
timeStep = 0.1; ithpower = 3;
options = struct('timeStep',timeStep,'ithpower',ithpower);
val = computeField(tls,'Apower_abs_dt_fact',options);
val_true = A_abs^3 * timeStep^3 / factorial(3);
assert(compareMatrices(val,val_true,tol,"equal",true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
