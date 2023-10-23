function res = testLong_linearSys_reach_07_constInput(~)
% testLong_linearSys_reach_07_constInput - unit test to check if constant inputs c
%  (cf. x' = Ax + Bu + c) are handled correctly
% note: the simulation results may be not with absolute certainty correct,
%       but should nonetheless remain inside the computed reachable sets
%
% Syntax:
%    testLong_linearSys_reach_07_constInput(~)
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false whether reachable sets overapproximative
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       07-January-2019
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% System dynamics ---------------------------------------------------------

A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linSys = linearSys('linearSys',A,B,c);

% Parameters --------------------------------------------------------------

params.R0 = zonotope([1;0],diag([0.05,0.05]));      % initial set
params.tFinal = 0.2;                                % final time
params.U = zonotope(0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.05;
options.taylorTerms = 10;
options.zonotopeOrder = 20;


% Reachability Analysis ---------------------------------------------------

R = reach(linSys,params,options);


% Simulation --------------------------------------------------------------

% number of initial points
simOpt.points = 5;
simRes = simulateRandom(linSys,params,simOpt); 


% Visualization -----------------------------------------------------------

plotting = false;
if plotting
    figure; hold on; box on;
    plot(R);
    plot(simRes);
end


% Numerical check ---------------------------------------------------------

% for every run, the simulation results must be within R
for i=1:simOpt.points
    disp("Trajectory " + i + "/" + simOpt.points + "...");
    results = simRes(i).x{1};
    
    for p=1:size(results,1)
        % loop over all points of the simulated trajectory
        p_curr = results(p,:)';
        
        for iSet=1:length(R.timeInterval.set)
            % check every time-interval set
            if contains(R.timeInterval.set{iSet},p_curr)
                break
            end
            if iSet == length(R.timeInterval.set)
                % simulated point not in computed reachable sets
                throw(CORAerror('CORA:testFailed'));
            end
        end
        
    end
end


end

% ------------------------------ END OF CODE ------------------------------
