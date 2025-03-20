function res = test_nn_rl_agentRL_validateOptions()
% test_nn_rl_agentRL_validateOptions - unit test function for
%     RLagents/@RLagent: checks coorect options validation.
%
% Syntax:
%    res = test_nn_rl_agentRL_validateOptions()
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

% initialize weights and biases of the layers
actorInit = dlnetwork(actorLayers);
nnActor = neuralNetwork.convertDLToolboxNetwork(num2cell(actorInit.Layers),false);

% initialize weights and biases of the layers
criticInit = dlnetwork(criticLayers);
nnCritic = neuralNetwork.convertDLToolboxNetwork(num2cell(criticInit.Layers),false);

% Test wrong combination for set-based critic
options.rl.critic.nn.train.method = 'naive';
options.rl.actor.nn.train.method = 'point';

assertThrowsAs(@agentDDPG,'CORA:wrongFieldValue',nnActor,nnCritic,options);

% Invalid parameters.
options.rl.critic.nn.train.method = 'grad';
options.rl.actor.nn.train.method = 'point';

assertThrowsAs(@agentDDPG,'CORA:wrongFieldValue',nnActor,nnCritic,options);

% Invalid parameters.
options.rl.critic.nn.train.method = 'set';
options.rl.actor.nn.train.method = 'point';

assertThrowsAs(@agentDDPG,'CORA:wrongFieldValue',nnActor,nnCritic,options);

% Invalid parameters.
options.rl.critic.nn.train.method = 'set';
options.rl.actor.nn.train.method = 'naive';

assertThrowsAs(@agentDDPG,'CORA:wrongFieldValue',nnActor,nnCritic,options);

options.rl.critic.nn.train.method = 'set';
options.rl.actor.nn.train.method = 'grad';

assertThrowsAs(@agentDDPG,'CORA:wrongFieldValue',nnActor,nnCritic,options);

% Test wrong type of propagation etchnique
options.rl.critic.nn.train.method = 'set';
options.rl.actor.nn.train.method = 'set';
options.rl.actor.nn.poly_method = 'regression';

assertThrowsAs(@agentDDPG,'CORA:wrongFieldValue',nnActor,nnCritic,options);

% Test wrong value for noise 
options.rl.actor.nn.poly_method = 'bounds';
options.rl.noise = -1;

assertThrowsAs(@agentDDPG,'CORA:outOfDomain',nnActor,nnCritic,options);

% Test wrong value for gamma 
options.rl.noise = .1;
options.rl.gamma = 2;

assertThrowsAs(@agentDDPG,'CORA:outOfDomain',nnActor,nnCritic,options);

% Test wrong Exploration Noise 
options.rl.gamma = .1;
options.rl.expNoiseType = 'Laplace';
options.rl.expNoise = -0.1;

assertThrowsAs(@agentDDPG,'CORA:outOfDomain',nnActor,nnCritic,options);

% Test wrong perturbation noise dimension
options.rl.expNoiseType = 'OU';
options.rl.expNoise = 0.1;
options.rl.noise = eye(2);

assertThrowsAs(@agentDDPG,'',nnActor,nnCritic,options);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
