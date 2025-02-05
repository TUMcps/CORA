function res = test_nonlinearSys_approxDepReach()
% test_nonlinearSys_approxDepReach - unit test function for approximative 
% reachabiltiy analysis of nonlinear systems (required for controller 
% synthesis)
%
% Syntax:
%    res = test_nonlinearSys_approxDepReach()
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

% Authors:       Lukas Schäfer
% Written:       07-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% dynamical system
tank = nonlinearSys(@tank6Eq); 

% parameters for reachability analysis
dim_x=6;
params.tFinal = 1; %final time
params.R0 = polyZonotope([2; 4; 4; 2; 10; 4],0.2*eye(dim_x));
params.U = zonotope([0,0.005]);

% algorithm settings for reachability analysis - 
options.zonotopeOrder = 1000; % choose a large value to avoid order reduction of the time-point reachable set
options.taylorTerms = 5;
options.tensorOrder = 3;
options.intermediateOrder = 10;
options.errorOrder = 10;
options.reductionTechnique = 'girard';
options.alg = 'poly';
options.timeStep = params.tFinal; % only one time step


% compute over-approximation of the reachable set
options.approxDepOnly = false;
R = reach(tank,params,options);

% compute approximation of the reachable set
options.approxDepOnly = true;
Rapprox = reach(tank,params,options);

% -> since the same parameters are used and we chose a large zonotope
% order (i.e., we do not need to account for reduction errors), we expect
% the approximative time-point reachable set to be a subset of the
% over-approximative reachable set - the same should apply for the
% respective zonotope enclosures since the additional generators of the
% zonotope representing the over-approximative reachable set are stored in 
% the independent generator matrix
assert(contains(zonotope(R.timePoint.set{end}),zonotope(Rapprox.timePoint.set{end}),'approx:st',1e-6));

% combine results
res = true;

end

% ------------------------------ END OF CODE ------------------------------
