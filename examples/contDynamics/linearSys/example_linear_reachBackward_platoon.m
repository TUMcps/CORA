function res = example_linear_reachBackward_platoon
% example_linear_reachBackward_platoon - example for backward reachability
%    analysis
%
% Syntax:
%    res = example_linear_reachBackward_platoon
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] F. Gruber and M. Althoff, "Computing Safe Sets of Linear
%        Sampled-Data Systems", IEEE Control Systems Letters 5 (2), 2021.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       20-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% target set at time 0
params.tStart = 0;
% solve backward until...
params.tFinal = 2;

% time step size and type of reachability (minimal, maximal)
options.timeStep = 0.02;
options.verbose = true;

% target set: goal set (EA)
% - relative distance should be above safe distance, but not too far
% - relative velocity should be within [-1.5,1.5]m/s (~5km/h)
% - acceleration should be within [-1,1]m/s^2 (~3.5km/h per second)
H_single_EA = [1 0 0; -1 0 0; -1 -2 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
d_single_EA = [20; 0; 0; 1.5; 1.5; 1; 1];

% target set: set to avoid (AE)
% - relative distance should not be below safe distance or crash
% - relative velocity should not be above 3m/s (~10km/h)
% - acceleration should not be above 1m/s^2 (~3.5km/h per second)
H_single_AE = [1 0 0; -1 0 0; -1 -2 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
d_single_AE = [0; 20; 7; 10; -3; 5; -1];

nrCarsLoop = [1,5,10];
R_EA_tp = cell(numel(nrCarsLoop),1);
R_EA_ti = cell(numel(nrCarsLoop),1);
R_AE_tp = cell(numel(nrCarsLoop),1);
R_AE_ti = cell(numel(nrCarsLoop),1);

for i=1:numel(nrCarsLoop)
    % number of cars
    nrCars = nrCarsLoop(i);

    % init platoon system
    [A,B,E] = platoon(nrCars);
    sys = linearSys(A,B,[],[],[],[],E);
    
    % set of controllable inputs:
    lower_bound = -5 * ones(nrCars,1);  % maximum deceleration: 5m/s^2
    upper_bound = 1 * ones(nrCars,1);   % maximum acceleration: 1m/s^2
    params.U = zonotope(interval(lower_bound,upper_bound));
    
    % disturbance set: acceleration of leading vehicle, see [1]
    % maximum acceleration and deceleration: +- 0.5 m/s^2
    params.W = zonotope(0,0.5);

    % target set for EA backward reachability
    H = kron(eye(nrCars),H_single_EA);
    d = repmat(d_single_EA,nrCars,1);
    params.R0 = polytope(H,d);
    
    % backward analyses
    options.linAlg = 'inner:EA:timepoint';
    R_EA_tp{i} = reachBackward(sys,params,options);
    options.linAlg = 'inner:EA:timeinterval';
    R_EA_ti{i} = reachBackward(sys,params,options);

    % target set for AE backward reachability
    H = kron(eye(nrCars),H_single_AE);
    d = repmat(d_single_AE,nrCars,1);
    params.R0 = polytope(H,d);

    options.linAlg = 'outer:AE:timepoint';
    R_AE_tp{i} = reachBackward(sys,params,options);
    options.linAlg = 'outer:AE:timeinterval';
    R_AE_ti{i} = reachBackward(sys,params,options);

end


% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
