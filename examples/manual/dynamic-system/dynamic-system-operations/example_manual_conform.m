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
params.tFinal = 5; % final time
params.R0conf = zonotope(zeros(2,1),eye(2)); % initial set
params.V = zonotope([0,1]); % sensor noise set
params.W = [-6; 1]*zonotope([0,1]); % disturbance set
y = [0.79; 5.00; 4.35; 1.86; -0.11; -1.13]; % measurement vector
delta = [0.1; 0.2; 0.1; 0; -0.1; -0.2]; % deviation of measurement vector
params.testSuite{1} = testCase(y, zeros(6,2), [1,1], 1); % test case 1
params.testSuite{2} = testCase(y + delta, zeros(6,2), [1,1], 1); % test case 2
params.testSuite{3} = testCase(y - delta, zeros(6,2), [1,1], 1); % test case 3

%% Algorithmic Settings
options.norm = 'interval'; % interval norm
options.timeStep = 1;

%% System Dynamics
reactor = linearSysDT('reactor',[0 -0.5; 1 1], 1, zeros(2,1), [-2 1], options.timeStep); 

% conformance synthesis
options.confAlg = 'dyn';
[params_conform, ~, ~, unifiedOutputs] = conform(reactor,params,options);

% plot --------------------------------------------------------------------

% Compute reachable set using obtained parameters
options.zonotopeOrder = inf;
params_conform.R0 = params_conform.R0conf;
R = reach(reactor, params_conform, options);

figure; hold on; box on;
useCORAcolors("CORA:contDynamics")
plotOverTime(R,1,'LineWidth',2);
updateColorIndex; % no initial set
% plot unified test cases
for iStep = 1:length(unifiedOutputs)
    plot(iStep-1, unifiedOutputs{iStep}(:,1),'k','Marker','.','LineStyle', 'none');
end

% ------------------------------ END OF CODE ------------------------------
