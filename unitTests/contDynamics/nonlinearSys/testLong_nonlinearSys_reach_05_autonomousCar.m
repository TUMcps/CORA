function res = testLong_nonlinearSys_reach_05_autonomousCar()
% testLong_nonlinearSys_reach_05_autonomousCar - unit_test_function of 
%    nonlinear reachability analysis for following a reference trajectory
%    Checks the solution of an autonomous car following a reference
%    trajectory; It is checked whether the reachable set is enclosed
%    in the initial set after a certain amount of time.
%
% Syntax:
%    res = testLong_nonlinearSys_reach_05_autonomousCar()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       10-September-2015
% Last update:   12-August-2016
%                23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

dim_x = 8;
params.tFinal=0.1; %final time
params.R0 = zonotope([[0; 0; 0; 22; 0 ; 0; -2.1854; 0],...
    0.05*diag([1, 1, 1, 1, 1, 1, 1, 1])]); %initial state for reachability analysis
params.u = uTRansVec4CASreach();
params.u = params.u(:,1:10);
params.U = zonotope([0*params.u(:,1), 0.05*diag([ones(5,1);zeros(21,1)])]);


% Reachability Settings ---------------------------------------------------

options.timeStep=0.01;
options.taylorTerms=5;
options.zonotopeOrder=200;
options.maxError = ones(dim_x,1); % for comparison reasons
options.alg = 'lin';
options.tensorOrder = 2;


% System Dynamics ---------------------------------------------------------

vehicle = nonlinearSys(@vmodel_A_bicycle_linear_controlled);


% Reachability Analysis --------------------------------------------------- 

R = reach(vehicle, params, options);


% Numerical Evaluation ----------------------------------------------------

%enclose result by interval
IH = interval(R.timeInterval.set{end});

%saved result
IH_saved = interval( ...
    [1.9113165972834527; -0.1763919418894342; -0.0628054820847352; 21.7144766641770168; -0.1081823785867376; -0.2077292813633349; -2.5375737787154540; -0.0418129308050226], ...
    [2.2486917066149412; 0.1775907076510275; 0.0638901147646950; 21.8628568496963389; 0.1329959142044251; 0.2499135613777886; -1.9819546869205054; 0.0646318759557788]);

%final result
assert(isequal(IH,IH_saved,1e-8));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
