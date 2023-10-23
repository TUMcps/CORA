function res = example_hybrid_reach_05_roomHeating()
% example_hybrid_reach_05_roomHeating - example for reachability of a
%    hybrid automaton considering the room heating benchmark described in 
%    Sec. 2.3 in [1] with two rooms
%
% Syntax:
%    example_hybrid_reach_05_roomHeating
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%   [1] A. Fehnker and F. Ivancic. "Benchmarks for Hybrid Systems 
%       Verification", HSCC 2004

% Authors:       Niklas Kochdumper
% Written:       26-June-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

HA = roomHeating();


% Parameter ---------------------------------------------------------------

params.tFinal = 12;
params.startLoc = 1;
params.R0 = zonotope([20.5;20.5],diag([0.1,0.1]));  
params.U = zonotope(4,0.01);


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.taylorTerms = 5; 
options.zonotopeOrder = 10; 
options.timeStep = 0.005; 
 
% settings for hybrid systems
options.enclose = {'box','pca'}; 
options.guardIntersect = 'zonoGirard';


% Simulation --------------------------------------------------------------

simRes = simulateRandom(HA,params);


% Reachability Analysis ---------------------------------------------------

R = reach(HA,params,options);


% Visualization -----------------------------------------------------------

% temperature room k over time
for k=1:2
    figure; box on; hold on;

    % plot reachable set
    plotOverTime(R,k);

    % plot simulation
    plotOverTime(simRes,k);

    % labels, title, and limitss
    xlabel('Time')
    ylabel('Temperature');
    title(['Room ' num2str(k)]);
    xlim([0,params.tFinal]);
end

res = true;

% ------------------------------ END OF CODE ------------------------------
