function res = test_linearSys_reach_10_zeroDyn()
% test_linearSys_reach_10_zeroDyn - unitTest to check if the trivial dynamics
% \dot{x} = 0 is correctly handled. If the inverse of the system matrix A=0
% would be required, this unit test would not pass.
%
% Syntax:
%    res = test_linearSys_reach_10_zeroDyn
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       30-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

dim_x = 5;
params.tFinal = 3; 
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = zonotope(0);


% Reachability Settings ---------------------------------------------------

options.timeStep      = 0.04;
options.taylorTerms   = 3;
options.zonotopeOrder = 100;


% System Dynamics ---------------------------------------------------------

A = zeros(dim_x);
B = zeros(dim_x,1);
zeroDynSys = linearSys('zeroDynSys',A,B);


% Reachability Analysis ---------------------------------------------------

options.linAlg = 'standard';
R_standard = reach(zeroDynSys, params, options);

options.linAlg = 'wrapping-free';
R_wrappingfree = reach(zeroDynSys, params, options);

options.linAlg = 'fromStart';
R_fromStart = reach(zeroDynSys, params, options);


% Is final set the same as initial set? -----------------------------------

% set accuracy to check if sets are matching
accuracy = 1e-8;

% standard approach
assert(isequal(R_standard.timePoint.set{end},params.R0,accuracy));

% wrapping-free approach
assert(isequal(R_wrappingfree.timePoint.set{end},params.R0,accuracy));

% fromStart approach
assert(isequal(R_fromStart.timePoint.set{end},params.R0,accuracy));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
