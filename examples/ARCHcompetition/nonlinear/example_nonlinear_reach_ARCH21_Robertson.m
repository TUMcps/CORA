function example_nonlinear_reach_ARCH21_Robertson()
% example_nonlinear_reach_ARCH21_Robertson - Robertson chemical reaction
%    example from the ARCH21 competition
%
% Syntax:  
%    example_nonlinear_reach_ARCH21_Robertson
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean  
% 
% Author:       Mark Wetzlinger
% Written:      17-May-2021
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

colors = {'b',[0,0.6,0],'r'};

% Parameters --------------------------------------------------------------

tFinal = 40;
params.tFinal = tFinal; 
params.U = zonotope(0);
R0 = zonotope([1;0;0]);

dim_x = dim(R0);
dim_u = dim(params.U);

% Reachability Settings ---------------------------------------------------

% Case 1: beta = 1e2, gamma = 1e3
options{1}.timeStep = 0.0025;
options{1}.taylorTerms = 6;
options{1}.zonotopeOrder = 30;
options{1}.intermediateOrder = 5;
options{1}.errorOrder = 5;

options{1}.alg = 'lin';
options{1}.tensorOrder = 3;

options{1}.verbose = true;

% Case 1: beta = 1e3, gamma = 1e5
options{2}.timeStep = 0.002;
options{2}.taylorTerms = 10;
options{2}.zonotopeOrder = 30;
options{2}.intermediateOrder = 5;
options{2}.errorOrder = 5;

options{2}.alg = 'lin';
options{2}.tensorOrder = 3;

% options{2}.alg = 'lin-adaptive';
options{2}.verbose = true;

% Case 1: beta = 1e3, gamma = 1e7
options{3}.timeStep = 0.00025;
options{3}.taylorTerms = 10;
options{3}.zonotopeOrder = 200;
options{3}.intermediateOrder = 5;
options{3}.errorOrder = 10;

options{3}.alg = 'lin';
options{3}.tensorOrder = 3;

% options{3}.alg = 'lin-adaptive';
options{3}.verbose = true;

% System Dynamics ---------------------------------------------------------

sys{1} = nonlinearSys(@Robertson_case1,dim_x,dim_u);
sys{2} = nonlinearSys(@Robertson_case2,dim_x,dim_u);
sys{3} = nonlinearSys(@Robertson_case3,dim_x,dim_u);

noSys = length(sys);

% Simulations -------------------------------------------------------------

% simOpt.points = 1;
% simOpt.fracVert = 1;
% simOpt.fracInpVert = 1;
% simOpt.inpChanges = 1;
% 
% for i=1:noSys
%     simRes{i} = simulateRandom(sys{i},params,simOpt);
% end
% 
% % plot over time
% figure;
% for d=1:dim_x
%     subplot(1,dim_x,d); hold on; box on;
%     xlabel("Time");
%     ylabel("x_" + d);
%     for i=1:noSys
%         han{i} = plotOverTime(simRes{i},d,'Color',colors{i});
%     end
%     legend('Case 1','Case 2','Case 3');
% end

% Reachability Analysis ---------------------------------------------------

R = cell(noSys,1);
compTime = zeros(noSys,1);
widthSum = cell(3,1);
volR = cell(3,1);

% split time horizon into two parts
tsplit = [5; 5; 2];
timeSteps = [0.0025, 0.025; 0.002, 0.005; 0.0002, 0.001];

currSys = 2;

for i=1:noSys %currSys:currSys %2:noSys
    
    params.R0 = R0;
    params.tFinal = tsplit(i);
    options{i}.timeStep = timeSteps(i,1);
    tic;
    R1{i} = reach(sys{i}, params, options{i});

    params.R0 = R1{i}.timePoint.set{end};
    params.tFinal = tFinal - params.tFinal;
    options{i}.timeStep = timeSteps(i,2);
    R2{i} = reach(sys{i}, params, options{i});
    % computation time
    compTime(i) = toc;

    % no time splitting
%     tic;
%     R{i} = reach(sys{i}, params, options{i});
%     % computation time
%     compTime(i) = toc;
    
    % append R2 to R1
    R{i} = append(R1{i},R2{i});
    
    % compute sum of concentrations of all time-interval sets
	noSets = length(R{i}.timeInterval.set);

    widthSum{i} = cell(noSets,1);
    for iSet = 1:noSets
        % width (interval) of each set
        widthSum{i}{iSet,1} = sum(interval(R{i}.timeInterval.set{iSet}));
    end
    
%     volR{i} = zeros(noSets,1);
%     for j=1:length(R{i}.timeInterval.set)
%         volR{i}(j,1) = volume(interval(R{i}.timeInterval.set{j}));
%     end
    
end



% Visualization -----------------------------------------------------------

% figure;
% projDims = {[1,2],[2,3]};
% for i=currSys:currSys
%     for p=1:2
%         subplot(1,2,p); hold on; box on;
%         plot(R{i},projDims{p},'FaceColor',colors{i},'EdgeColor',colors{i},...
%             'Filled',true);
%         xlabel("x_" + projDims{p}(1));
%         ylabel("x_" + projDims{p}(2));
%     end
% end


figure; hold on; box on;

% plot sum of concentrations over time
for i=currSys:currSys %1:noSys
    Rtime = R{i}.timeInterval.time;

    for j=1:length(widthSum{i})

        % get intervals
        intX = widthSum{i}{j};
        intT = Rtime{j};

        int = cartProd(intT,intX);

        % plot interval
        han{i} = plot(int,[1,2],'FaceColor',colors{i},...
            'EdgeColor','none','Filled',true);
    end
    
end
    
% label plot
xlabel('t');
ylabel('s');
% axis([0,params.tFinal,0.999,1.001]);
legend('Case 1', 'Case 2', 'Case 3');

% close all;

%------------- END OF CODE --------------