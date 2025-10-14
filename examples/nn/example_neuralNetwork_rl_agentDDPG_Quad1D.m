function completed = example_neuralNetwork_rl_agentDDPG_Quad1D()
% example_neuralNetwork_rl_agentDDPG_Quad1D - Example for set-based
%   reinforcement learning. This example script trains a point-based
%   and set-based [1]. The results are only evaluated for one
%   random seed. In the example we use a shortened evaluation. To fully
%   train the agents set 'totalEpisodes' to 2000 and initialOps to 'uniform'.
%
% Syntax:
%    res = example_neuralNetwork_rl_agentDDPG_Quad1D()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% References:
%   [1] Wendl, M. et al. Training Verifiably Robust Agents Using Set-Based
%       Reinforcement Learning, 2024.
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Manuel Wendl, Tobias Ladner
% Written:       27-August-2024
% Last update:   12-October-2024 (TL, shortened example, marked long runs as benchmark)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------
m = 0.05;
g = 9.81;
% open-loop system
f = @(x,u) [x(2);(u(1)+1)/(2*m)-g];
sys = nonlinearSys(f);

% Settings ----------------------------------------------------------------

totalEpisodes = 300;

% RL settings
options.rl.gamma = .99;
options.rl.tau = .005;
options.rl.expNoise = .2;
options.rl.expNoiseType = 'OU';
options.rl.expDecay = 1;
options.rl.batchsize = 64;
options.rl.buffersize = 1e6;
options.rl.noise = .1;
options.rl.printFreq = 10;
options.rl.visRate = 10;

% Actor Settings
options.rl.actor.nn.use_approx_error = true;
options.rl.actor.nn.poly_method = 'bounds';
options.rl.actor.nn.train.optim = nnAdamOptimizer(1e-4,.9,.999,1e-8,0);
options.rl.actor.nn.train.method = 'point';
options.rl.actor.nn.train.eta = .1;
options.rl.actor.nn.train.omega = 0;
options.rl.actor.nn.train.zonotope_weight_update = 'sum';
options.rl.actor.nn.train.backprop = true;
options.rl.actor.nn.train.advOps.numSamples = 200;
options.rl.actor.nn.train.advOps.alpha = 4;
options.rl.actor.nn.train.advOps.beta = 4;

% Critic Settings
options.rl.critic.nn.use_approx_error = true;
options.rl.critic.nn.poly_method = 'bounds';
options.rl.critic.nn.train.optim = nnAdamOptimizer(1e-3,.9,.999,1e-8,1e-2);
options.rl.critic.nn.train.method = 'point';
options.rl.critic.nn.train.eta = .01;
options.rl.critic.nn.train.zonotope_weight_update = 'sum';
options.rl.critic.nn.train.backprop = true;

% Env settings
options.rl.env.x0 = interval([-4;0],[4;0]);
options.rl.env.initialOps = 'symmetric';
options.rl.env.dt = .1;
options.rl.env.timeStep = .01;
options.rl.env.maxSteps = 30;
options.rl.env.collisioCheckBool = true;
options.rl.env.evalMode = 'point';
options.rl.env.solver = 'ODE45';
options.rl.env.reach.alg = 'lin';
options.rl.env.reach.zonotopeOrder = 200;

% Build Environment -------------------------------------------------------

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);

% Build NNs ---------------------------------------------------------------

% Actor NN 
actorLayers = [
    featureInputLayer(2)
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(1)
    tanhLayer
    ];

% Critic NN 
criticLayers = [
    featureInputLayer(3)
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(1)
    ];

seed = 1;
rng(seed,"twister"); % Set rng
% initialize weights and biases of the layers
actorInit = dlnetwork(actorLayers);
nnActor = neuralNetwork.convertDLToolboxNetwork(num2cell(actorInit.Layers),false);

% initialize weights and biases of the layers
criticInit = dlnetwork(criticLayers);
nnCritic = neuralNetwork.convertDLToolboxNetwork(num2cell(criticInit.Layers),false);

% Train networks ----------------------------------------------------------

% Point-Based Actor and Point-Based Critic ---
options.rl.actor.nn.train.method = 'point';
options.rl.critic.nn.train.method = 'point';

rng(seed,"twister"); % Set rng
DDPG_PAPC = agentDDPG(nnActor,nnCritic,options);
DDPG_PAPC = DDPG_PAPC.train(env,totalEpisodes);

% Set-Based Actor and Set-Based Critic ---
options.rl.actor.nn.train.method = 'set';
options.rl.actor.nn.train.omega = .5;
options.rl.critic.nn.train.method = 'set';

rng(seed,"twister"); % Set rng
DDPG_SASC = agentDDPG(nnActor,nnCritic,options);
DDPG_SASC = DDPG_SASC.train(env,totalEpisodes);

% Evaluation --------------------------------------------------------------

% Parameters ---

params.R0 = zonotope(interval([3.8;0],[4.2;0]));
params.tFinal = 5;

% Reachability settings ---

options_eval.timeStep = 0.1;
options_eval.alg = 'lin';
options_eval.tensorOrder = 2;
options_eval.taylorTerms = 4;
options_eval.zonotopeOrder = 100;

% Options for NN evaluation ---

options_eval.nn = struct();
options_eval.nn.poly_method = "singh";

% compute reachable sets ---

% PAPC
sys_PAPC = neurNetContrSys(sys,DDPG_PAPC.actor.nn,options.rl.env.dt);
R_PAPC = reach(sys_PAPC,params,options_eval);

% SASC
sys_SASC = neurNetContrSys(sys,DDPG_SASC.actor.nn,options.rl.env.dt);
R_SASC = reach(sys_SASC,params,options_eval);

% Visualization -----------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:contDynamics",2)
plotOverTime(R_PAPC,1,'DisplayName','Reachable set PA-PC')
plotOverTime(R_SASC,1,'DisplayName','Reachable set SA-SC')
legend(); xlabel('Time'); ylabel('x')

completed = true;

end


% Auxiliary functions -----------------------------------------------------

function reward = aux_rewardFun_Quadrocopter1D(x)
% reward function
if isnumeric(x)
    reward = -(abs(x(end,1))+0.01*abs(x(end,2)));
else
    w = [1,0.01,0];
    reward = -max(abs(supportFunc(x.timePoint.set{end},w,'lower')),abs(supportFunc(x.timePoint.set{end},w)));
end
end

function collisionBool = aux_collisionCheck_Quadrocopter1D(x)
% collision check
if isa(x,'numeric')
    if all(abs(x(:,1:2))<0.05,"all")
        collisionBool = true;
    else
        collisionBool = false;
    end
elseif isa(x,'reachSet')
    if norm(x.timePoint.set{1}.c)<0.05
        collisionBool = true;
    else
        collisionBool = false;
    end
end
end

% ------------------------------ END OF CODE ------------------------------
