function res = testnn_rl_agentRL_train_point()
% testnn_rl_agentRL_train_point - unit test function for 
%     RLagent/train: train a point-based RLagent with CORA and 
%       check if reward is greater than a given threshold after 
%       75 learning itterations.
%
% Syntax:
%    res = testnn_rl_agentRL_train_point()
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% System Dynamics ---------------------------------------------------------
m = 0.05;
g = 9.81;
% open-loop system
f = @(x,u) [x(2);(u(1)+1)/(2*m)-g];
sys = nonlinearSys(f);

% Actor NN ----------------------------------------------------------------
actorLayers = [
    featureInputLayer(2)
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(1)
    tanhLayer
    ];

% Critic NN ---------------------------------------------------------------
criticLayers = [
    featureInputLayer(3)
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(1)
    ];

% Settings ----------------------------------------------------------------
totalEpisodes = 100;

% RL settings
options.rl.gamma = .99;
options.rl.tau = .005;
options.rl.printFreq = 1;
options.rl.expNoise = 0.1;

% Env settings
options.env.dt = .1;
options.env.timeStep = .01;
options.env.maxSteps = 30;
options.env.x0 = interval([-4;0],[4;0]);
options.env.initialOps = 'uniform';

% Build Environment -------------------------------------------------------

env = ctrlEnvironment(sys,@aux_rewardFun_Quadrocopter1D,@aux_collisionCheck_Quadrocopter1D,options);

% Build NNs ---------------------------------------------------------------

rng(1,"twister"); % Set rng
% initialize weights and biases of the layers
actorInit = dlnetwork(actorLayers);
nnActor = neuralNetwork.convertDLToolboxNetwork(num2cell(actorInit.Layers),false);

% initialize weights and biases of the layers
criticInit = dlnetwork(criticLayers);
nnCritic = neuralNetwork.convertDLToolboxNetwork(num2cell(criticInit.Layers),false);

% Point based Agent -------------------------------------------------------
DDPG1 = agentDDPG(nnActor,nnCritic,options);
DDPG1 = DDPG1.train(env,totalEpisodes,true);

%% Evaluate Learning Histories

r = movmean(DDPG1.learnHistory.reward,10);

% Implementation reward > threshold of -10
assert(mean(r(75:end)) > -10);

end


% Auxiliary functions -----------------------------------------------------

function reward = aux_rewardFun_Quadrocopter1D(x)
% reward function of quadrocopter
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
