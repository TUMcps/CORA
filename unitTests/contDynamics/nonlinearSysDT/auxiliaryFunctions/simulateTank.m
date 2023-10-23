function simulateTank(options)
% simulateTank - simulates a tank system to disturbances and 
% sensor noise
%
% Syntax:
%    simulateTank(options)
%
% Inputs:
%    options - options struct
%
% Outputs:
%    sys - system struct
%    options - options struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Matthias Althoff
% Written:       25-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
savepath = [CORAROOT filesep 'unitTests' filesep 'contDynamics' ...
    filesep 'nonlinearSysDT' filesep 'models'];

%% Settings

options.timeStep = 0.5;
options.zonotopeOrder = 100;
options.reductionTechnique = 'pca';
options.tensorOrder = 2;
options.errorOrder = 1;

%% Input trajectory 
% We use the change of rate of the input to sepcify trajectories

% input patterns
segmentDuration = 10; %[s]
segmentSteps = ceil(segmentDuration/options.timeStep);
uDiff = 5e-3*[...
    1 2 3 2  1  0 -1 -3 -3 0 -1 -2 -1; ...
    1 1 1 0  0  0 -1 -1 -1 0 -1 -1  1;
    2 1 1 1  0 -1 -1  0 -1 0  0  0  1];
    %1 2 1 0 -1 -1 -1 -1  0 1  1  0  0]; 
% input
u = 0.1*ones(length(uDiff(:,1))); 
for i = 1:length(uDiff(1,:))
    for iInter = 1:segmentSteps
        u(:,end+1) = u(:,end) + options.timeStep*uDiff(:,i);
    end
end

% select measurement matrix
C = zeros(3,6);
C(1,2) = 1; % second tank measured
C(2,4) = 1; % fourth tank measured
C(3,5) = 1; % fifth tank measured

% dimension, inputs, and measurements
dim = length(C(1,:));
measurements = length(C(:,1));
inputs = length(u(:,1));

%% System Dynamics 

% helper function
fun = @(x,u) tank6EqDT_inflow4(x,u,options.timeStep);

tank = nonlinearSysDT('tankSystem',fun,options.timeStep,dim,inputs,C); % initialize tank system


%% Parameters
W = zonotope([zeros(dim,1), 0.001*eye(dim)]); % disturbance set
V = zonotope([zeros(measurements,1), 0.2*eye(measurements)]); % sensor noise set

params.tStart = 0; %start time --> remove after fixing options checks
params.tFinal = length(u(1,:))*options.timeStep; %final time
params.R0 = zonotope(20*ones(dim,1),4*eye(dim)); %initial set
params.V = V; % sensor noise set
params.W = W; % disturbance set
params.u = u; %input transition


%% Simulation Settings

options.points = 1;
options.p_conf = 0.999; % probability that sample of normal distribution within specified set
options.type = 'gaussian';

% simulate result assuming Gaussian distributions
simRes = simulateRandom(tank, params, options);

%% obtain output values
for i=1:length(simRes.t{1})
    % create measurement noise
    v = randPoint(params.V,'gaussian',options.p_conf);
    % obtain output value
    params.y(:,i) = C*simRes.x{1}(i,:)' + v;
end

save([savepath '/' 'tankModel_nonlin_dim6_new'], 'tank', 'params', 'options', 'simRes');

% ------------------------------ END OF CODE ------------------------------
