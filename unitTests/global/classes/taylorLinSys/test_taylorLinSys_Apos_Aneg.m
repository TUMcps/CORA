function res = test_taylorLinSys_Apos_Aneg
% test_taylorLinSys_Apos_Aneg - unit test for helper class taylorLinSys
%
% Syntax:
%    res = test_taylorLinSys_Apos_Aneg
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

% compare to directly computed solution
options = struct('ithpower',1);
val = computeField(tls,'Apos',options);
val_true = A; val_true(val_true < 0) = 0;
assert(compareMatrices(val,val_true,tol,"equal",true));

% compute higher power where there are gaps in between
options = struct('ithpower',5);
val = computeField(tls,'Apos',options);
val_true = A^5; val_true(val_true < 0) = 0;
assert(compareMatrices(val,val_true,tol,"equal",true));


% start with power higher than 1
options = struct('ithpower',3);
val = computeField(tls,'Aneg',options);
val_true = A^3; val_true(val_true > 0) = 0;
assert(compareMatrices(val,val_true,tol,"equal",true));

% compute next-higher power
options = struct('ithpower',4);
val = computeField(tls,'Aneg',options);
val_true = A^4; val_true(val_true > 0) = 0;
assert(compareMatrices(val,val_true,tol,"equal",true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
