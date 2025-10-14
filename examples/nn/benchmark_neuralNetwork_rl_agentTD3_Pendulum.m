function completed = benchmark_neuralNetwork_rl_agentTD3_Pendulum()
% benchmark_neuralNetwork_rl_agentTD3_Pendulum - Example for set-based
%   reinforcement learning. This example script trains a point-based,
%   naive-adversarial, grad-adversarial, SA-PC and SA-SC agents with
%   omega = 0 and omega = 0.5 [1]. The results are only evaluated for one
%   random seed. To fully train the agents set 'totalEpisodes' to 2000. 
%
% Syntax:
%    res = benchmark_neuralNetwork_rl_agentTD3_Pendulum()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
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

% Authors:       Manuel Wendl
% Written:       27-August-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------
m = 1;
l = 1;
g = 9.81;
% open-loop system
f = @(x,u) [x(2);g/l*sin(x(1))+1/(m*l^2)*15*u(1)];
sys = nonlinearSys(f);

% Settings ----------------------------------------------------------------
totalEpisodes = 300;

% RL settings
options.rl.gamma = .99;
options.rl.tau = .005;
options.rl.expNoise = .1;
options.rl.expNoiseTarget = .2;
options.rl.expNoiseType = 'gaussian';
options.rl.expDecay = 1;
options.rl.batchsize = 64;
options.rl.buffersize = 1e6;
options.rl.noise = .1;
options.rl.printFreq = 25;
options.rl.visRate = 25;

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
options.rl.env.x0 = interval([-pi;0],[pi;0]);
options.rl.env.initialOps = 'uniform';
options.rl.env.dt = .1;
options.rl.env.timeStep = .01;
options.rl.env.maxSteps = 30;
options.rl.env.collisioCheckBool = true;
options.rl.env.evalMode = 'point';
options.rl.env.solver = 'ODE45';
options.rl.env.reach.alg = 'lin';
options.rl.env.reach.zonotopeOrder = 200;

% Allocate Agent Array ----------------------------------------------------
seed = 1;
numAgents = 6;
agents = cell(1,numAgents);

% Build Environment -------------------------------------------------------

env = ctrlEnvironment(sys,@aux_rewardFun_Pendulum,@aux_collisionCheck_Pendulum,options);

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

rng(seed,"twister"); % Set rng
% initialize weights and biases of the layers
actorInit = dlnetwork(actorLayers);
nnActor = neuralNetwork.convertDLToolboxNetwork(num2cell(actorInit.Layers),false);

% initialize weights and biases of the layers
criticInit1 = dlnetwork(criticLayers);
nnCritic1 = neuralNetwork.convertDLToolboxNetwork(num2cell(criticInit1.Layers),false);

criticInit2 = dlnetwork(criticLayers);
nnCritic2 = neuralNetwork.convertDLToolboxNetwork(num2cell(criticInit2.Layers),false);

% Point-Based Actor -------------------------------------------------------
options.rl.actor.nn.train.method = 'point';
options.rl.critic.nn.train.method = 'point';

rng(seed,"twister"); % Set rng
TD3_1 = agentTD3(nnActor,nnCritic1,nnCritic2,options);
TD3_1 = TD3_1.train(env,totalEpisodes);
agents{1} = TD3_1;

% Point-Based Actor naive adversarial attack ------------------------------
options.rl.actor.nn.train.method = 'naive';
options.rl.critic.nn.train.method = 'point';

rng(seed,"twister"); % Set rng
TD3_2 = agentTD3(nnActor,nnCritic1,nnCritic2,options);
TD3_2 = TD3_2.train(env,totalEpisodes);
agents{2} = TD3_2;

% Point-Based Actor grad adversarial attack -------------------------------
options.rl.actor.nn.train.method = 'grad';
options.rl.critic.nn.train.method = 'point';

rng(seed,"twister"); % Set rng
TD3_3 = agentTD3(nnActor,nnCritic1,nnCritic2,options);
TD3_3 = TD3_3.train(env,totalEpisodes);
agents{3} = TD3_3;

% Set-Based Actor and Point-Based Critic ----------------------------------
options.rl.actor.nn.train.method = 'set';
options.rl.critic.nn.train.method = 'point';

rng(seed,"twister"); % Set rng
TD3_4 = agentTD3(nnActor,nnCritic1,nnCritic2,options);
TD3_4 = TD3_4.train(env,totalEpisodes);
agents{4} = TD3_4;

% Set-Based Actor and Set-Based Critic ------------------------------------
options.rl.actor.nn.train.method = 'set';
options.rl.actor.nn.train.omega = 0;
options.rl.critic.nn.train.method = 'set';

rng(seed,"twister"); % Set rng
TD3_5 = agentTD3(nnActor,nnCritic1,nnCritic2,options);
TD3_5 = TD3_5.train(env,totalEpisodes);
agents{5} = TD3_5;

% Set-Based Actor and Set-Based Critic ------------------------------------
options.rl.actor.nn.train.method = 'set';
options.rl.actor.nn.train.omega = .5;
options.rl.critic.nn.train.method = 'set';

rng(seed,"twister"); % Set rng
TD3_6 = agentTD3(nnActor,nnCritic1,nnCritic2,options);
TD3_6 = TD3_6.train(env,totalEpisodes);
agents{6} = TD3_6;

[~,Reward,RewardAdvNaive,RewardAdvGrad] = compareAgents(agents,env,0:0.01:0.1,{'z','dz','a'},'sup');

%% Evaluate Learning Histories

meanLbReward = cellfun(@mean,Reward);
meanNaiveReward = cellfun(@mean,RewardAdvNaive);
meanGradReward = cellfun(@mean,RewardAdvGrad);

disp('Mean LB Rewards:')
disp(meanLbReward)

disp('Mean Naive Rewards:')
disp(meanNaiveReward)

disp('Mean Grad Rewards:')
disp(meanGradReward)

completed = true;

end


% Auxiliary functions -----------------------------------------------------

function reward = aux_rewardFun_Pendulum(x)
% reward function
if isnumeric(x)
    reward = -(abs(x(end,1))+0.01*abs(x(end,2)));
else
    w = [1,0.01,0];
    reward = -max(abs(supportFunc(x.timePoint.set{end},w,'lower')),abs(supportFunc(x.timePoint.set{end},w)));
end
end

function collisionBool = aux_collisionCheck_Pendulum(x)
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
