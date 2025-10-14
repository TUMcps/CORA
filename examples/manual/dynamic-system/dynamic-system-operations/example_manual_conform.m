function example_manual_conform()
% example_manual_conform - example from the manual demonstrating the 
% conformance operations as defined in the manual
%
% Syntax:
%   example_manual_conform()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Matthias Althoff
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Parameters
dt = 1; 
params.tFinal = 5; % final time
params.R0 = zonotope(zeros(2,1),eye(2)); % initial set
params.V = zonotope([0,1]); % sensor noise set
params.W = [-6; 1]*zonotope([0,1]); % disturbance set
y = [0.79 5.00 4.35 1.86 -0.11 -1.13]; %measurement vector
delta = [0.1 0.2 0.1 0 -0.1 -0.2]; %deviation of measurement vector
params.testSuite = [trajectory(zeros(2,6), [1;1], y, [], dt); %test case 1
    trajectory(zeros(2,6), [1;1], y + delta, [], dt); %test case 2
    trajectory(zeros(2,6), [1;1], y - delta, [], dt)];

%% Algorithmic Settings
options.norm = 'interval'; % interval norm
options.timeStep = 1;

%% System Dynamics
A = [0 -0.5; 1 1];
B = [1; 1];
C = [-2 1];
E = eye(dim(params.W));
F = eye(dim(params.V));
reactor = linearSysDT('reactor',A,B,[],C,[],[],E,F,options.timeStep); 

% conformance synthesis
[params_conform, results] = conform(reactor,params,options);
unifiedOutputs = results.unifiedOutputs;

% plot --------------------------------------------------------------------

% Compute reachable set using obtained parameters
options.zonotopeOrder = inf;
R = reach(reactor, params_conform, options);

figure; hold on; box on;
useCORAcolors("CORA:contDynamics")
plotOverTime(R,1,'LineWidth',2);
updateColorIndex; % no initial set
% plot unified test cases
for iStep = 1:size(unifiedOutputs,2)
    plot(iStep-1, squeeze(unifiedOutputs(1,iStep,:)),'k','Marker','.','LineStyle', 'none');
end

% ------------------------------ END OF CODE ------------------------------
