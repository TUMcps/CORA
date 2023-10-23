function simulateTank_30(options)
% simulateTank_30 - simulates a tank system to disturbances and 
% sensor noise
%
% Syntax:
%    simulateTank_30(options)
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
% Written:       28-March-2021
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
% uDiff = 5e-3*[...
%     1 2 3 2  1  0 -1 -3 -3 0 -1 -2 -1; ...
%     1 1 1 0  0  0 -1 -1 -1 0 -1 -1  1; ...
%     2 1 1 1  0 -1 -1  0 -1 0  0  0  1; ...
%     1 2 3 2  1  0 -1 -3 -3 0 -1 -2 -1; ...
%     1 1 1 0  0  0 -1 -1 -1 0 -1 -1  1; ...
%     2 1 1 1  0 -1 -1  0 -1 0  0  0  1; ...
%     1 2 3 2  1  0 -1 -3 -3 0 -1 -2 -1; ...
%     1 1 1 0  0  0 -1 -1 -1 0 -1 -1  1; ...
%     2 1 1 1  0 -1 -1  0 -1 0  0  0  1; ...
%     1 2 3 2  1  0 -1 -3 -3 0 -1 -2 -1; ...
%     1 1 1 0  0  0 -1 -1 -1 0 -1 -1  1; ...
%     2 1 1 1  0 -1 -1  0 -1 0  0  0  1; ...
%     1 2 3 2  1  0 -1 -3 -3 0 -1 -2 -1; ...
%     1 1 1 0  0  0 -1 -1 -1 0 -1 -1  1; ...
%     2 1 1 1  0 -1 -1  0 -1 0  0  0  1];
uDiff = 5e-3*[...
    1 2 3 2  1  0 -1 -3 -3 0 -1; ...
    1 1 1 0  0  0 -1 -1 -1 0 -1; ...
    2 1 1 1  0 -1 -1  0 -1 0  0; ...
    1 2 3 2  1  0 -1 -3 -3 0 -1; ...
    1 1 1 0  0  0 -1 -1 -1 0 -1; ...
    2 1 1 1  0 -1 -1  0 -1 0  0; ...
    1 2 3 2  1  0 -1 -3 -3 0 -1; ...
    1 1 1 0  0  0 -1 -1 -1 0 -1; ...
    2 1 1 1  0 -1 -1  0 -1 0  0; ...
    1 2 3 2  1  0 -1 -3 -3 0 -1; ...
    1 1 1 0  0  0 -1 -1 -1 0 -1; ...
    2 1 1 1  0 -1 -1  0 -1 0  0; ...
    1 2 3 2  1  0 -1 -3 -3 0 -1; ...
    1 1 1 0  0  0 -1 -1 -1 0 -1; ...
    2 1 1 1  0 -1 -1  0 -1 0  0];
% input
u = 0.1*ones(length(uDiff(:,1))); 
for i = 1:length(uDiff(1,:))
    for iInter = 1:segmentSteps
        u(:,end+1) = u(:,end) + options.timeStep*uDiff(:,i);
    end
end

% select measurement matrix
C = zeros(21,30);
C(1,2) = 1; % second tank measured
C(2,4) = 1; % fourth tank measured
C(3,5) = 1; % fifth tank measured
C(4,7) = 1; % 7th tank measured
C(5,8) = 1; % 8th tank measured
C(6,10) = 1; % 10th tank measured
C(7,11) = 1; % 11th tank measured
C(8,13) = 1; % 13th tank measured
C(9,14) = 1; % 14th tank measured
C(10,16) = 1; % 16th tank measured
C(11,17) = 1; % 17th tank measured
C(12,19) = 1; % 19th tank measured
C(13,20) = 1; % 20th tank measured
C(14,21) = 1; % 21th tank measured
C(15,22) = 1; % 22th tank measured
C(16,23) = 1; % 23th tank measured
C(17,25) = 1; % 25th tank measured
C(18,26) = 1; % 26th tank measured
C(19,27) = 1; % 27th tank measured
C(20,28) = 1; % 28th tank measured
C(21,29) = 1; % 29th tank measured

% dimension, inputs, and measurements
n = length(C(1,:));
measurements = length(C(:,1));
inputs = length(u(:,1));

out_fun = @(x,u) C*x(1:n);
outputs = size(C,1);

%% System Dynamics 

% helper function
fun = @(x,u) tank30EqDT_inflow15(x,u,options.timeStep);

tank = nonlinearSysDT('tankSystem_30',fun,options.timeStep,n,inputs,out_fun,outputs); % initialize tank system


%% Parameters
W = zonotope([zeros(n,1), 0.001*eye(n)]); % disturbance set
V = zonotope([zeros(measurements,1), 0.2*eye(measurements)]); % sensor noise set

params.tFinal = length(u(1,:))*options.timeStep; %final time
params.R0 = zonotope(20*ones(n,1),4*eye(n)); %initial set
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

% change input set due to prototypical implementation
params.U = zeros(inputs,1); % input set

save([savepath '/' 'tankModel_nonlin_dim30_new'], 'tank', 'params', 'options', 'simRes');

% ------------------------------ END OF CODE ------------------------------
