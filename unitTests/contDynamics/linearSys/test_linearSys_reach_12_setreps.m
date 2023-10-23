function res = test_linearSys_reach_12_setreps()
% test_linearSys_reach_12_setreps - unit test function of linear
%    reachability analysis with uncertain inputs for different set
%    representations for the initial set and uncertain inputs
%
% Syntax:
%    res = test_linearSys_reach_12_setreps()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       17-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

U = interval(-0.01,0.02);
params.tFinal = 0.1;


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05; % time step size for reachable set computation
options.taylorTerms = 4; % number of taylor terms in exponential matrix
options.zonotopeOrder = 10; % zonotope order
options.linAlg = 'standard';


% System Dynamics ---------------------------------------------------------

A = [-1 -4; 4 -1];
B = [1; 0.5];
sys = linearSys('twoDimSys',A,B);


% Reachability Analysis ---------------------------------------------------

% zonotopes
params.R0 = zonotope([10;5],[0.5 0.1; -0.4 0.2]);
params.U = zonotope(U);
R = reach(sys, params, options);

% constrained zonotope
params.R0 = conZonotope([10;5],[0.2 0.1 -0.3; -0.1 0.3 0.2],[0.3 -0.5 0.2],0);
params.U = conZonotope(U);
R = reach(sys, params, options);

% polyZonotope
params.R0 = polyZonotope([10;5],[0.2 0.1; -0.1 0.3],[0.2;0.1],[2 0; 0 1]);
params.U = polyZonotope(U);
R = reach(sys, params, options);

% ellipsoid
params.R0 = ellipsoid(zonotope([10;5],[0.5 0.1; -0.4 0.2]));
params.U = ellipsoid(U);
R = reach(sys, params, options);

% capsule
params.R0 = capsule([10;5],[0.2;-0.1],0.1);
params.U = capsule(U);
R = reach(sys, params, options);

% polytope: no convex hull function...
% params.R0 = [10;5] + polytope([1 1; 0 -1; 1 -1],ones(3,1));
% params.U = polytope(U);
% R = reach(sys, params, options);

% all executed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
