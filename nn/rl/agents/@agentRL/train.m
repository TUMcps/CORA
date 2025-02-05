function obj = train(obj,env,varargin)
% train - train agentRL with initial state
%
% Syntax:
%   obj = train(obj,env,varargin)
%
% Inputs:
%   obj - agentRL
%   env - control environment
%   episodes - number of train episodes
%   verbose - boolean to print training log in terminal
% 
% Outputs:
%   obj - trained DDPGagent
%   learnHistory - learninig history 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: DDPGagent, actor, critic

% Authors:       Manuel Wendl
% Written:       24-October-2023
% Last update:   18-August-2024 (extract environment from RlAgent, new options structure)
%                19-September-2024 (TL, renamed to train)
%                12-October-2024 (TL, speed up set-based exp. gathering)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parse function arguments
narginchk(2,4)
[episodes, verbose] = setDefaultValues({2000, true},varargin);

inputArgsCheck({ ...
    {env, 'att', 'ctrlEnvironment'}; ...
    {episodes, 'att', 'numeric'}; ...
    {verbose, 'att', {'logical','numeric'}};
    })

table = aux_printInfo(verbose,obj.options,env.options.rl.env,episodes);
learnHistory = aux_initLearningHistory(episodes);

% Get inputDataClass
% To speed up computations and reduce gpu memory, we only use single 
% precision.
inputDataClass = single(1);
% Check if a gpu is used during training.
useGpu = obj.options.rl.actor.nn.train.use_gpu;
if useGpu
    % Training data is also moved to gpu.
    inputDataClass = gpuArray(inputDataClass);
end


% Allocate generators for initial perturbance set.
if numel(obj.options.rl.noise) == 1
    % If pertubation radius is given
    idMat = eye(obj.options.rl.actor.nn.train.num_init_gens,'like',inputDataClass);
    noiseBatchG = cast(repmat(idMat*obj.options.rl.noise,1,1,obj.options.rl.batchsize),'like',inputDataClass);
else
    % if custom pertubation generators are given check
    % dimensions
    if size(obj.options.rl.noise,1) ~= obj.environment.ctrlDynamics.nrOfDims
        throw(CORAerror("CORA:dimensionMismatch",obj.options.rl.noise,obj.environment.ctrlDynamics.nrOfDims));
    end
    noiseBatchG = cast(repmat(obj.options.rl.noise,1,1,obj.options.rl.batchsize),'like',inputDataClass);
    obj.options.rl.noise = mean(sum(abs(obj.options.rl.noise),2));
end

% Reset buffer from previous learning
obj.buffer = obj.buffer.resetBuffer();

% initialise exploration factor for exploration decay
explorationFactor = 1;

% start timer
tic

% Begin Training Loop -----------------------------------------------------
for episode = 1:episodes

    % Get initial observation of each episode
    [env, observation] = env.reset();
    observation = cast(observation,'like',inputDataClass);
    isDone = false;

    % Update exploration Decay
    explorationFactor = aux_updateExplorationFactor(explorationFactor,obj.options);

    % Initial exploration noise 
    e = zeros(obj.actor.nn.neurons_out,1,'like',inputDataClass);

    % Initialize update counter for reward normalization in logging
    numberOfUpdates = 1;
    
    while ~isDone
        
        % Gather Experience -----------------------------------------------

        % evaluate actor 
        if strcmp(obj.options.rl.critic.nn.train.method,'point')
            if any(strcmp(obj.options.rl.actor.nn.train.method,{'point','set'}))
                action = obj.actor.nn.evaluate_(observation,obj.options.rl.actor,obj.actor.idxLayer);
            else
                % Perform adversarial attack if not point or set-based
                advObs = obj.computeAdversarialAttack(observation);
                action = obj.actor.nn.evaluate_(advObs,obj.options.rl.actor,obj.actor.idxLayer);
            end
        elseif strcmp(obj.options.rl.critic.nn.train.method,'set')
            [action,actionG] = obj.actor.nn.evaluateZonotopeBatch_(observation,idMat*obj.options.rl.noise,obj.options.rl.actor,obj.actor.idxLayer);
            zAction = zonotope(action,actionG);
        else
            throw(CORAerror("CORA:notDefined",['There are no other ' ...
                'training methods than point and set-based training ' ...
                'used for the critic.']))
        end
        
        % add exploration noise
        e = aux_getExplorationNoise(e,obj.options);
        action = action + explorationFactor*e;
        % limit action
        if ~isa(obj.actor.nn.layers{end},'nnLinearLayer')
            action = max(min(action,1),-1);
        end

        % evaluate environment
        [env, nextObservation, reward, isDone, visual] = env.step(action);
        
        % cast environemnt outputs
        nextObservation = cast(nextObservation,'like',inputDataClass);
        isDone = cast(isDone,'like',inputDataClass);
        reward = cast(reward,'like',inputDataClass);

        % store transition in replay buffer
        if strcmp(obj.options.rl.critic.nn.train.method,'set')
            data = {observation,zAction,reward,nextObservation,isDone};
        else
            data = {observation,action,reward,nextObservation,isDone};
        end
        obj.buffer = obj.buffer.fillBuffer(data,obj.options);

        % log Q0 and reward
        if env.stepNum == 2
            Q0 = obj.critic.nn.evaluate_([observation;action],obj.options.rl.critic,obj.critic.idxLayer);
            learnHistory.Q0(episode) = Q0;
        end
        
        % Set current state to nextState
        observation = nextObservation;
        
        % Accumulate reward
        learnHistory.reward(episode) = learnHistory.reward(episode) + reward;

        % log visual 
        if mod(episode,obj.options.rl.visRate) == 0
            obj.buffer = obj.buffer.storeVisualData(visual,episode);
        end

        % Train Networks --------------------------------------------------

        % Checkif enough experience was gathered
        if obj.buffer.currentIndx > obj.options.rl.batchsize || numberOfUpdates > 1
            % get random batch
            randBatch = obj.buffer.getRandomBatch(obj.options);
            
            if isa(obj,'TD3agent') % In case of TD3 add target exp. noise
                 randBatch{4} = randBatch{4} + randn(size(randBatch{4}),'like',randBatch{4})*obj.options.rl.expNoiseTarget;
            end

            % perform network update step
            [obj, learnHistory] = obj.trainNetworkStep(randBatch,noiseBatchG,learnHistory,episode);

            numberOfUpdates = numberOfUpdates + 1;
        end
    end
    
    learnHistory = aux_normalizeLossHistory(learnHistory,episode,numberOfUpdates);
    
    time = toc;
    
    if 0 == mod(episode,obj.options.rl.printFreq)
        aux_updateInfo(table,learnHistory,time,episode,verbose)
    end

    if episode > obj.options.rl.earlyStop
        if std(learnHistory.reward(episode-obj.options.rl.earlyStop:episode)) < abs(0.01*mean(learnHistory.reward(episode-obj.options.rl.earlyStop:episode)))
            if verbose
                fprintf(['Stopped early! The reward for the last %i ' ...
                    'episodes did not increase.\n'], ...
                        obj.options.rl.earlyStop);
            end
            learnHistory.criticLoss.center = learnHistory.criticLoss.center(1:episode);
            learnHistory.criticLoss.vol = learnHistory.criticLoss.vol(1:episode);
            learnHistory.actorLoss.center = learnHistory.actorLoss.center(1:episode);
            learnHistory.actorLoss.vol = learnHistory.actorLoss.vol(1:episode);
            learnHistory.reward = learnHistory.reward(1:episode);
            learnHistory.Q0 = learnHistory.Q0(1:episode);
            obj.learnHistory = learnHistory;
            break;
        end
    end
    
end

if verbose
    table.printFooter();    
end

obj.learnHistory = learnHistory;

% (Potentially) gather weights of networks from GPU 
obj.gatherNetworks();

% Delete all gradients
obj.deleteAllGradients();

end


% Auxiliary functions -----------------------------------------------------

function table = aux_printInfo(verbose,options,envOps,episodes)
if verbose
    % init table
    table = CORAtableParameters('Reinforcement Learning Parameters');
    table.printHeader();

    % print standard rl parameters
    table.printContentRow('Standard-RL Parameter')
    table.printContentRow('Discount Factor (gamma)', options.rl.gamma,'.2e',2)
    table.printContentRow('Target update factor (tau)', options.rl.tau,'.2e',2)
    table.printContentRow('Exploration noise type', options.rl.expNoiseType,'s',2)
    table.printContentRow('Exploration noise', options.rl.expNoise,'.2e',2)
    table.printContentRow('Batch size', options.rl.batchsize,'i',2)
    table.printContentRow('Buffer size', options.rl.buffersize,'i',2)
    table.printContentRow('Episodes', episodes,'i',2)
    
    % actor parameters
    if strcmp(options.rl.actor.nn.train.method,'set') && strcmp(options.rl.critic.nn.train.method,'set')
        % print set-based actor (& critic) parameter
        table.printContentRow('Set-based Actor Training')
        table.printContentRow('Noise',options.rl.noise,'.2e',2)
        table.printContentRow('Actor weighting factor (eta)',options.rl.actor.nn.train.eta,'.2e',2)
        table.printContentRow('Prior weighting factor (omega)',options.rl.actor.nn.train.omega,'.2e',2)
    elseif strcmp(options.rl.actor.nn.train.method,'set')
        % print set-based actor parameter
        table.printContentRow('Set-based Actor Training')
        table.printContentRow('Noise',options.rl.noise,'.2e',2)
        table.printContentRow('Actor weighting factor (eta)',options.rl.actor.nn.train.eta,'.2e',2)
    elseif strcmp(options.rl.actor.nn.train.method,'point')
        % print point-based actor parameter
        table.printContentRow('Point-based Actor Training')
    else
        % print adversarial actor parameter
        table.printContentRow('Adversarial Actor Training')
        table.printContentRow('Method',options.rl.actor.nn.train.method,'s',2)
    end
    table.printContentRow('Optimizer',options.rl.actor.nn.train.optim.print(),'s',2)
    
    % critic parameters
    if strcmp(options.rl.critic.nn.train.method,'set')
        table.printContentRow('Set-based Critic Training')
        table.printContentRow('Critic weighting factor (eta)',options.rl.critic.nn.train.eta,'.2e',2)
    elseif strcmp(options.rl.critic.nn.train.method,'point')
        table.printContentRow('Point-based Critic Training')
    end
    table.printContentRow('Optimizer',options.rl.critic.nn.train.optim.print(),'s',2)

    % environment parameters
    table.printContentRow('Environment Parameters')
    table.printContentRow('Control time step (dt)',envOps.dt,'.2e',2)
    table.printContentRow('Simulation time step',envOps.timeStep,'.2e',2)
    table.printContentRow('Max control steps/epoch',envOps.maxSteps,'.2e',2)
    table.printContentRow('Point solver',envOps.solver,'s',2)
    table.printContentRow('Set alg',envOps.reach.alg,'s',2)
    table.printContentRow('Tensor order',envOps.reach.tensorOrder,'i',2)
    table.printContentRow('Taylor terms',envOps.reach.taylorTerms,'i',2)
    table.printContentRow('Max zonotope order',envOps.reach.zonotopeOrder,'i',2)

    % finish parameter table
    table.printFooter();

    % start training table
    hvalues = {'Epoch','Training Time','Actor Loss','Critic Loss','Reward','Q0'};
    formats = {'d','s','.2f','.2f','.2f','.2f'};
    colWidths = [0,0,0,0,10,10];
    table = CORAtable('double',hvalues,formats,'ColumnWidths',colWidths);
    table.printHeader();
else
    table = [];
end

end

function aux_updateInfo(table,learnHistory,trainTime,episode,verbose)
if verbose
    trainTimeVec = [0 0 0 0 0 trainTime];
    cvalues = {episode,datetime(trainTimeVec,'Format','HH:mm:ss'), ...
        learnHistory.actorLoss.center(episode), learnHistory.criticLoss.center(episode), ...
        learnHistory.reward(episode), learnHistory.Q0(episode) ...
    };
    table.printContentRow(cvalues);
end
end

function learnHistory = aux_initLearningHistory(episodes)
% Initialise learn history
learnHistory.criticLoss.center = zeros(episodes,1);
learnHistory.criticLoss.vol = zeros(episodes,1);
learnHistory.actorLoss.center = zeros(episodes,1);
learnHistory.actorLoss.vol = zeros(episodes,1);
learnHistory.reward = zeros(episodes,1);
learnHistory.Q0 = zeros(episodes,1);
end

function explorationFactor = aux_updateExplorationFactor(explorationFactor,options)
if options.rl.expDecayFactor > 0
    if options.rl.expDecayFactor > 1
        throw(CORAerror("CORA:wrongValue","options.rl.decayFactor","Possible Decay factors for exponential decay are [0,1] and for linear decay are < 0."))
    end
    explorationFactor = explorationFactor * options.rl.expDecayFactor;
else
    explorationFactor = max(0,explorationFactor + options.rl.expDecayFactor);
end
end

function e = aux_getExplorationNoise(e,options)
% add random gaussian exploration noise
if strcmp(options.rl.expNoiseType,'gaussian')
    e = randn(size(e),'like',e)*options.rl.expNoise;
elseif strcmp(options.rl.expNoiseType,'OU')
% add OU Noise
    e = e - 0.15*e*.01 + options.rl.expNoise*sqrt(.01)*randn(size(e),'like',e);
else
    throw(CORAerror("CORA:specialError",'This noise type is not known.'))
end
end

function learnHistory = aux_normalizeLossHistory(learnHistory,episode,numberOfUpdates)
learnHistory.criticLoss.center(episode) = learnHistory.criticLoss.center(episode)/numberOfUpdates;
learnHistory.criticLoss.vol(episode) = learnHistory.criticLoss.vol(episode)/numberOfUpdates;
learnHistory.actorLoss.center(episode) = learnHistory.actorLoss.center(episode)/numberOfUpdates;
learnHistory.actorLoss.vol(episode) = learnHistory.actorLoss.vol(episode)/numberOfUpdates;
end

% ------------------------------ END OF CODE ------------------------------
