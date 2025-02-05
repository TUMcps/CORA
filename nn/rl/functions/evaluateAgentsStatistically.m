function evaluateAgentsStatistically(agents,env,perturbation,initialOps)
% evaluateAgentsStatistically - statistical evaluation of runs with different
% random seed
%
% Syntax:
%   evaluateAgentsStatistically(agents,env,perturbation,initialOps)
%
% Inputs:
%   agents - cell Array of trained agents of size (seed,different RLagents)
%        All RLagents in one row have to have the same rng seed!
%   env - control environment
%   pertubation - input pertubation radii (array)
%   initialOps - inintial observation option for evaluation
%
% Outputs:
%   None
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: RLagent, ctrlEnvironment

% Authors:       Manuel Wendl
% Written:       26-November-2023
% Last update:   18-August-2024 (extract environment from RLagent, new options structure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% This script changes all interpreters from tex to latex.
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

% This script forces box on
set(groot,'DefaultAxesBox','on')

numberOfSeeds = size(agents,1);
numberofAgents = size(agents,2);

% Learning History Allocation
LHReward = zeros(numberOfSeeds,numberofAgents,length(agents{1,1}.learnHistory.reward));
LHQ0 = zeros(numberOfSeeds,numberofAgents,length(agents{1,1}.learnHistory.Q0));

% Evaluation Measures Allocation
Reward = zeros(numberOfSeeds,numberofAgents,length(perturbation));
RewardAdvNaive = zeros(numberOfSeeds,numberofAgents,length(perturbation));
RewardAdvGrad = zeros(numberOfSeeds,numberofAgents,length(perturbation));

for i = 1:numberOfSeeds
    for j = 1:numberofAgents
        reward = zeros(size(perturbation));
        rewardAdvNaive = zeros(size(perturbation));
        rewardAdvGrad = zeros(size(perturbation));

        for k = 1:length(perturbation)
            [~,reward(k),rewardAdvNaive(k),rewardAdvGrad(k)] = agents{i,j}.benchmark(env,perturbation(k),initialOps);
        end

        LHReward(i,j,:) = agents{i,j}.learnHistory.reward;
        LHQ0(i,j,:) = agents{i,j}.learnHistory.Q0;

        Reward(i,j,:) = reward;
        RewardAdvNaive(i,j,:) = rewardAdvNaive;
        RewardAdvGrad(i,j,:) = rewardAdvGrad;
    end
end

% Specifiy plot colors and labels for each agent
cmap = zeros(numberofAgents,3);
labels = cell(1,numberofAgents);

for j = 1:numberofAgents
    if strcmp(agents{1,j}.options.rl.actor.nn.train.method,'point') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'point')
        cmap(j,:) = [0.3,0.3,0.3];
        labels{j} = 'PA-PC';
    elseif strcmp(agents{1,j}.options.rl.actor.nn.train.method,'rand') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'point')
        cmap(j,:) = [14,255,0]./255;
        labels{j} = 'Random';
        % gradient-based
    elseif strcmp(agents{1,j}.options.rl.actor.nn.train.method,'extreme') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'point')
        cmap(j,:) = [31,198,0]./255;
        labels{j} = 'Extreme';
    elseif strcmp(agents{1,j}.options.rl.actor.nn.train.method,'naive') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'point')
        cmap(j,:) = [8,144,0]./255;
        labels{j} = 'Naive';
    elseif strcmp(agents{1,j}.options.rl.actor.nn.train.method,'grad') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'point')
        cmap(j,:) = [10,93,0]./255;
        labels{j} = 'Grad';
        % set-based 
    elseif strcmp(agents{1,j}.options.rl.actor.nn.train.method,'set') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'point')
        cmap(j,:) = [0.69020,0.82350,1.00000];
        labels{j} = 'SA-PC';
    elseif strcmp(agents{1,j}.options.rl.actor.nn.train.method,'set') && strcmp(agents{1,j}.options.rl.critic.nn.train.method,'set')
        if agents{1,j}.options.rl.actor.nn.train.omega == 0
            cmap(j,:) = [0.15300,0.41960,0.72740];
            labels{j} = 'SA-SC-$\omega=0$';
        else
            cmap(j,:) = [0.03530,0.25100,0.45490];
            labels{j} = 'SA-SC-$\omega=0.5$';
        end
    end
end

% prepare rewards
meanLHReward = squeeze(mean(LHReward,1));
confLHReward = tinv(.95,numberOfSeeds-1)*squeeze(std(LHReward,0,1))/sqrt(numberOfSeeds);

meanLHQ0 = squeeze(mean(LHQ0,1));

meanReward = squeeze(mean(Reward,1));
confReward = tinv(.95,numberOfSeeds-1)*squeeze(std(Reward,0,1))/sqrt(numberOfSeeds);

meanRewardAdvNaive = squeeze(mean(RewardAdvNaive,1));
confRewardAdvNaive = tinv(.95,numberOfSeeds-1)*squeeze(std(RewardAdvNaive,0,1))/sqrt(numberOfSeeds);

meanRewardAdvGrad = squeeze(mean(RewardAdvGrad,1));
confRewardAdvGrad = tinv(.95,numberOfSeeds-1)*squeeze(std(RewardAdvGrad,0,1))/sqrt(numberOfSeeds);

% visualization ---

% plot LH reward
fig1 = figure();
hold on
for i = 1:numberofAgents
    fill([1:size(meanLHReward,2), size(meanLHReward,2):-1:1],[meanLHReward(i,:)+confLHReward(i,:), meanLHReward(i,end:-1:1)-confLHReward(i,end:-1:1)],cmap(i,:),'FaceAlpha',0.2,'LineStyle','none','DisplayName',[labels{i},' - ','Learn History Reward 0.95 confidence']);
    plot(1:size(meanLHReward,2),meanLHReward(i,:),'LineStyle','-','Color',cmap(i,:),'DisplayName',[labels{i},' - ','Mean Learn History Reward'])
    plot(1:size(meanLHReward,2),meanLHQ0(i,:),'LineStyle','--','Color',cmap(i,:),'DisplayName',[labels{i},' - ','Mean Learn History Q0'])
end
xlabel('Episode')
ylabel('Reward during Training')
legend
savefig(fig1,'statisticLearnHistory.fig');

% plot rewards
fig2 = figure();
% plot mean reward
subplot(1,3,1)
hold on
for i = 1:numberofAgents
    fill([perturbation, perturbation(end:-1:1)],[meanReward(i,:)+confReward(i,:), meanReward(i,end:-1:1)-confReward(i,end:-1:1)],cmap(i,:),'FaceAlpha',0.2,'LineStyle','none','DisplayName',[labels{i},' - ',' 0.95 confidence']);
    plot(perturbation,meanReward(i,:),'Color',cmap(i,:),'DisplayName',labels{i})
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('LB')
legend
% plot adv naive reward
subplot(1,3,2)
hold on
for i = 1:numberofAgents
    fill([perturbation, perturbation(end:-1:1)],[meanRewardAdvNaive(i,:)+confRewardAdvNaive(i,:), meanRewardAdvNaive(i,end:-1:1)-confRewardAdvNaive(i,end:-1:1)],cmap(i,:),'FaceAlpha',0.2,'LineStyle','none','DisplayName',[labels{i},' - ',' 0.95 confidence']);
    plot(perturbation,meanRewardAdvNaive(i,:),'Color',cmap(i,:),'DisplayName',labels{i})
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('Naive')
% plot adv grad reward
subplot(1,3,3)
hold on
for i = 1:numberofAgents
    fill([perturbation, perturbation(end:-1:1)],[meanRewardAdvGrad(i,:)+confRewardAdvGrad(i,:), meanRewardAdvGrad(i,end:-1:1)-confRewardAdvGrad(i,end:-1:1)],cmap(i,:),'FaceAlpha',0.2,'LineStyle','none');
    plot(perturbation,meanRewardAdvGrad(i,:),'Color',cmap(i,:))
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
title('Grad')
savefig(fig2,'statisticEvaluation.fig');

% plot lower bound
fig3 = figure();
hold on
for i = 1:numberofAgents
    fill([perturbation, perturbation(end:-1:1)],[meanReward(i,:)+confReward(i,:), meanReward(i,end:-1:1)-confReward(i,end:-1:1)],cmap(i,:),'FaceAlpha',0.2,'LineStyle','none','DisplayName',[labels{i},' - ',' 0.95 confidence']);
    plot(perturbation,meanReward(i,:),'Color',cmap(i,:),'DisplayName',labels{i})
end
xlabel('$\epsilon$')
ylabel('$\underline{V}_\mu(s(t_0))$')
legend
savefig(fig3,'statisticEvaluationLB.fig');


% ------------------------------ END OF CODE ------------------------------
