function [reachSets,Reward,RewardAdvNaive,RewardAdvGrad] = compareAgents(agents,env,pertubation,stateLabels,initialOps)
% compareAgents - comparison of point based and set based agents after
% training.
%
% Syntax:
%   [reachSets,Reward,RewardAdvNaive,RewardAdvGrad] = compareAgents(Agents,env,pertubation,stateLabels,initialOps)
%
% Inputs:
%   agents - cell array of trained agents
%   env - control environment
%   pertubation - input pertubation radii (array)
%   stateLabels - cell array with state labels 
%   initialOps - initialOps for evaluation
%
% Outputs:
%   None
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ctrlEnvironment

% Authors:       Manuel Wendl
% Written:       15-November-2023
% Last update:   18-August-2024 (extract environment from RLagent)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

numAgents = length(agents);

% Do not save figures if called by test functions
saveBool = true;
stack = dbstack('-completenames');
if numel(stack) > 1
    if contains(stack(2).name,'test_')||contains(stack(2).name,'example_')
        saveBool = false;
    end
end

% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

% This script forces box on
set(groot,'DefaultAxesBox','on')

linestyles = {'-','--','-.',':','-','--'};

% Specify plot colors and labels for each agent
cmap = zeros(numAgents,3);
labels = cell(1,numAgents);

for i = 1:numAgents
    if strcmp(agents{i}.options.rl.actor.nn.train.method,'point') && strcmp(agents{i}.options.rl.critic.nn.train.method,'point')
        cmap(i,:) = [0.3,0.3,0.3];
        labels{i} = 'PA-PC';
    elseif strcmp(agents{i}.options.rl.actor.nn.train.method,'rand') && strcmp(agents{i}.options.rl.critic.nn.train.method,'point')
        cmap(i,:) = [14,255,0]./255;
        labels{i} = 'Random';
        % gradient-based
    elseif strcmp(agents{i}.options.rl.actor.nn.train.method,'extreme') && strcmp(agents{i}.options.rl.critic.nn.train.method,'point')
        cmap(i,:) = [31,198,0]./255;
        labels{i} = 'Extreme';
    elseif strcmp(agents{i}.options.rl.actor.nn.train.method,'naive') && strcmp(agents{i}.options.rl.critic.nn.train.method,'point')
        cmap(i,:) = [8,144,0]./255;
        labels{i} = 'Naive';
    elseif strcmp(agents{i}.options.rl.actor.nn.train.method,'grad') && strcmp(agents{i}.options.rl.critic.nn.train.method,'point')
        cmap(i,:) = [10,93,0]./255;
        labels{i} = 'Grad';
        % set-based
    elseif strcmp(agents{i}.options.rl.actor.nn.train.method,'set') && strcmp(agents{i}.options.rl.critic.nn.train.method,'point')
        cmap(i,:) = [0.69020,0.82350,1.00000];
        labels{i} = 'SA-PC';
    elseif strcmp(agents{i}.options.rl.actor.nn.train.method,'set') && strcmp(agents{i}.options.rl.critic.nn.train.method,'set')
        if agents{i}.options.rl.actor.nn.train.omega == 0
            cmap(i,:) = [0.15300,0.41960,0.72740];
            labels{i} = 'SA-SC-$\omega=0$';
        else
            cmap(i,:) = [0.03530,0.25100,0.45490];
            labels{i} = 'SA-SC-$\omega=0.5$';
        end
    end
end

% Plot trajectories 
fig_trajectories = figure();
for i = 1:4
    for j = 1:numAgents
        ind = i*floor(length(agents{j}.buffer.visualisationData.episodeNum)/4);
        % Create a subplot
        subplot(4,numAgents,numAgents*(i-1)+j)
        hold on
        for k = 1:length(stateLabels)
            plot(agents{j}.buffer.visualisationData.reachSet{ind}.t,agents{j}.buffer.visualisationData.reachSet{ind}.x(:,k),'Color',cmap(j,:),'LineStyle',linestyles{k},'DisplayName',[labels{j},': ',stateLabels{k}])
        end
        ylabel(strcat('Episode: ',num2str(agents{j}.buffer.visualisationData.episodeNum(ind))))
        xlabel('$t [s]$')
        if i == 1
            legend
        end
    end
end
if saveBool
    savefig(fig_trajectories,'Trajectories.fig');
end

% Plot Learning History 
fig_learingHistory = figure();
subplot(3,1,1) % reward
hold on
for i = 1:numAgents
    plot(agents{i}.learnHistory.reward,'Color',cmap(i,:),'LineStyle','-','DisplayName',[labels{i},': $\sum r_i$'])
    plot(agents{i}.learnHistory.Q0,'Color',cmap(i,:),'LineStyle','--','DisplayName',[labels{i},': $Q_0$'])
end
xlabel('Epsiode')
ylabel('$\sum r_i$/$Q_0$')
legend

subplot(3,1,2) % critic loss
hold on
for i = 1:numAgents
    yyaxis left
    plot(agents{i}.learnHistory.criticLoss.center,'Color',cmap(i,:),'LineStyle','-','Marker','none','DisplayName',[labels{i},': $\overline{\mathcal{L}_{Q,center}}$'])
    ylabel('$\overline{\mathcal{L}_{Q,center}}$')
    if strcmp(agents{i}.options.rl.critic.nn.train.method,'set')
        yyaxis right
        plot(agents{i}.learnHistory.criticLoss.vol,'Color',cmap(i,:),'LineStyle','--','DisplayName',[labels{i},': $\overline{\mathcal{L}_{Q,vol}}$'])
        ylabel('$\overline{\mathcal{L}_{Q,vol}}$')
    end
end
xlabel('Epsiode')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend

subplot(3,1,3) % actor loss
hold on
for i = 1:numAgents
    yyaxis left
    % Plot the loss
    plot(agents{i}.learnHistory.actorLoss.center,'Color',cmap(i,:),'LineStyle','-','Marker','none','DisplayName',[labels{i},': $\overline{\mathcal{L}_{a,center}}$'])
    ylabel('$\overline{\mathcal{L}_{a,center}}$')
    if strcmp(agents{i}.options.rl.actor.nn.train.method,'set')
        yyaxis right
        plot(agents{i}.learnHistory.actorLoss.vol,'Color',cmap(i,:),'LineStyle','--','DisplayName',[labels{i},': $\overline{\mathcal{L}_{a,vol}}$'])
        ylabel('$\overline{\mathcal{L}_{a,vol}}$')
    end
end
xlabel('Episode')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
legend

if saveBool
    savefig(fig_learingHistory,'LearningHistory.fig');
end

% Evaluate Volume Losses and Reachsets for differnet Pertubations
Reward = cell(1,numAgents);
RewardAdvNaive = cell(1,numAgents);
RewardAdvGrad = cell(1,numAgents);

for i = 1:numAgents
    Reward{i} = zeros(size(pertubation));
    RewardAdvNaive{i} = zeros(size(pertubation));
    RewardAdvGrad{i} = zeros(size(pertubation));
end

reachSets = cell(1,numAgents);

ind = floor(linspace(1,length(pertubation),4));

% plot reachable sets
fig_ReachSets = figure();
plotNum = 1;
for i = 1:length(pertubation)
    volSort = zeros(1,numAgents);
    for j = 1:numAgents
        % evaluate each agent
        [reachSet,reward,rewardAdvNaive,rewardAdvGrad] = agents{j}.benchmark(env,pertubation(i),initialOps);
        
        volSort(j) = sum(abs(reachSet(end).timePoint.set{end}.G),'all');
        Reward{j}(i) = reward;
        RewardAdvNaive{j}(i) = rewardAdvNaive;
        RewardAdvGrad{j}(i) = rewardAdvGrad;
        reachSets{j} = reachSet;
    end

    % plot in descending order to not cover each other
    [~,sortInd] = sort(volSort,'descend');
    if any(i == ind)
        for k = 1:length(stateLabels)
            subplot(4,length(stateLabels),plotNum+k-1)
            hold on
            for j = 1:numAgents
                try
                    plotOverTime(reachSets{sortInd(j)},k,'FaceColor',cmap(sortInd(j),:),'Unify',true,'DisplayName',labels{sortInd(j)})
                catch ME
                    fprintf('Unable to plot reachable set of agent %i: %s\n', j, ME.message);
                end
            end
            if plotNum == 1
                legend
            end
            xlabel('$t [s]$')
            ylabel(stateLabels{k})
            %ylim([-6,6])
            title(strcat('$\epsilon = $',num2str(pertubation(i))))
        end
        plotNum = plotNum + length(stateLabels);
    end
end
if saveBool
    savefig(fig_ReachSets,'ReachSets.fig')
end

% plot lower bound
fig_ComparisonLowerBound = figure();
hold on 
for i = 1:numAgents
    plot(pertubation,Reward{i},'Color',cmap(i,:),'DisplayName',labels{i})
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('LB')
legend
ylim([min(Reward{i}(Reward{i}>-1e8)),0])
if saveBool
    savefig(fig_ComparisonLowerBound,'ComparisonLowerBound.fig')
end

% plot lower bound attacks
fig_ComparisonLowerBoundAndAttacks = figure();
subplot(1,3,1) % reward
hold on 
minReward = 0;
for i = 1:numAgents
    plot(pertubation,Reward{i},'Color',cmap(i,:),'DisplayName',labels{i})
    minReward = min(minReward,min(Reward{i}(Reward{i}>-1e8)));
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('LB')
legend
ylim([minReward,0])

subplot(1,3,2) % reward adv naive
hold on 
minReward = 0;
for i = 1:numAgents
    plot(pertubation,RewardAdvNaive{i},'Color',cmap(i,:),'DisplayName',labels{i})
    minReward = min(minReward,min(Reward{i}(Reward{i}>-1e8)));
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('Adv. Naive Return')
legend
ylim([minReward,0])

subplot(1,3,3) % reward adv grad
hold on 
for i = 1:numAgents
    plot(pertubation,RewardAdvGrad{i},'Color',cmap(i,:),'DisplayName',labels{i})
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('Adv. Grad Return')
legend
if saveBool
    savefig(fig_ComparisonLowerBoundAndAttacks,'ComparisonLowerBoundAndAttacks.fig')
end
end

% ------------------------------ END OF CODE ------------------------------
