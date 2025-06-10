function res = testnn_rl_agentRL_train_all()
% testnn_rl_agentRL_train_all - unit test function for rlAgent/train: 
%   point-based, naive-adversarial, gradient-adversarial, set-based (SA-PC, SA-SC)
%   (only checks if code runs through)
%
% Syntax:
%    res = testnn_rl_agentRL_train_all()
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
% See also:

% Authors:       Manuel Wendl
% Written:       27-August-2024
% Last update:   08-May-2025 (TL, reduced train time, removed comparison, clean up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set seed for random initialization of networks
rng(4286746681);

% System Dynamics ---------------------------------------------------------
m = 0.05;
g = 9.81;
% open-loop system
f = @(x,u) [x(2);(u(1)+1)/(2*m)-g];
sys = nonlinearSys(f);

% Settings ----------------------------------------------------------------
totalEpisodes = 4;

% RL settings
options.rl.gamma = .99;
options.rl.tau = .005;
options.rl.expNoise = .2;
options.rl.expNoiseType = 'OU';
options.rl.expDecay = 1;
options.rl.batchsize = 64;
options.rl.buffersize = 1e6;
options.rl.noise = .1;
options.rl.printFreq = 1;
options.rl.visRate = 1;

% Actor Settings
options.rl.actor.nn.use_approx_error = true;
options.rl.actor.nn.poly_method = 'bounds';
options.rl.actor.nn.train.optim = nnAdamOptimizer(1e-4,.9,.999,1e-8,0);
options.rl.actor.nn.train.method = 'point';
options.rl.actor.nn.train.eta = .1;
options.rl.actor.nn.train.omega = 0;
options.rl.actor.nn.train.zonotope_weight_update = 'outer_product';
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
options.rl.critic.nn.train.zonotope_weight_update = 'outer_product';
options.rl.critic.nn.train.backprop = true;

% Env settings
options.rl.env.x0 = interval([-4;0],[4;0]);
options.rl.env.initialOps = 'uniform';
options.rl.env.dt = .1;
options.rl.env.timeStep = .1;
options.rl.env.maxSteps = 30;
options.rl.env.collisioCheckBool = true;
options.rl.env.evalMode = 'point';
options.rl.env.solver = 'ODE45';
options.rl.env.reach.alg = 'lin';
options.rl.env.reach.zonotopeOrder = 20;

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

rng(1,"twister"); % Set rng
% initialize weights and biases of the layers
actorInit = dlnetwork(actorLayers);
nnActor = neuralNetwork.convertDLToolboxNetwork(num2cell(actorInit.Layers),false);

% initialize weights and biases of the layers
criticInit = dlnetwork(criticLayers);
nnCritic = neuralNetwork.convertDLToolboxNetwork(num2cell(criticInit.Layers),false);

% Allocate Agent Array ----------------------------------------------------
seed = 1;
numAgents = 2;
agents = cell(1,numAgents);

% Point based Agent naive advers attack ---
options.rl.actor.nn.train.method = 'naive';
options.rl.critic.nn.train.method = 'point';

rng(seed,"twister"); % Set rng
DDPG1 = agentDDPG(nnActor,nnCritic,options);
DDPG1 = DDPG1.train(env,totalEpisodes);
agents{1} = DDPG1;

% Set based Agent and set based Critic ---
options.rl.actor.nn.train.method = 'set';
options.rl.actor.nn.train.omega = 0.5;
options.rl.critic.nn.train.method = 'set';

rng(seed,"twister"); % Set rng
DDPG2 = agentDDPG(nnActor,nnCritic,options);
DDPG2 = DDPG2.train(env,totalEpisodes);
agents{2} = DDPG2;

% evaluate ---
disp('Evaluate actors ...')
[~,Reward,RewardAdvNaive,RewardAdvGrad] = compareAgents(agents,env,0:0.1:0.2,{'z','dz','a'},'inf');

% clean up files
delete ComparisonLowerBound.fig; close;
delete ComparisonLowerBoundAndAttacks.fig; close;
delete LearningHistory.fig; close;
delete ReachSets.fig; close;
delete Trajectories.fig; close;

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function reward = aux_rewardFun_Quadrocopter1D(x)
% reward functions of quadrocopter
if isnumeric(x)
    reward = -(abs(x(end,1))+0.01*abs(x(end,2)));
else
    w = [1,0.01,0];
    reward = -max(abs(supportFunc(x.timePoint.set{end},w,'lower')),abs(supportFunc(x.timePoint.set{end},w)));
end
end

function collisionBool = aux_collisionCheck_Quadrocopter1D(x)
% collision check function
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
