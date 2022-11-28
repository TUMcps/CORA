function text = example_nonlinear_reach_ARCH22_Robertson()
% example_nonlinear_reach_ARCH22_Robertson - Robertson chemical reaction
%    example from the ARCH22 competition
%
% Syntax:  
%    example_nonlinear_reach_ARCH22_Robertson
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%

% Author:       Mark Wetzlinger
% Written:      17-May-2021
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

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

nrSys = length(sys);


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

R = cell(nrSys,1);
compTime = zeros(nrSys,1);
widthSum = cell(nrSys,1);
widthEnd = zeros(nrSys,1);

% split time horizon into two parts
tsplit = [5; 5; 2];
timeSteps = [0.0025, 0.025; 0.002, 0.005; 0.0002, 0.001];

for i=1:nrSys % 0
    
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
    widthEnd(i,1) = 2*rad(widthSum{i}{end});
    
end

% spare re-computation
% load ROBE22.mat

% for i=1:nrSys
%     for iSet = 1:length(R{i}.timeInterval.set)
%         widthSum{i}{iSet,1} = sum(interval(R{i}.timeInterval.set{iSet}));
%     end
% end



% Visualization -----------------------------------------------------------

% caution: plotting of R{3} takes a long time!

colors = {colorblind('b'),colorblind('r'),colorblind('y')};
figure; hold on; box on;

% init reachSet object for widthSum (faster plotting)
Rwidth = cell(nrSys,1);
for i=1:nrSys
    % time and time-interval solution
    timeInt.time = R{i}.timeInterval.time;
    timeInt.set = widthSum{i};
    Rwidth{i} = reachSet([],timeInt);

    % plot
    plotOverTime(Rwidth{i},1,'FaceColor',colors{i},'Unify',true);

    % flush queue
    pause(0.2);
end

% label plot
xlabel('t');
ylabel('s');
% axis([0,params.tFinal,0.999,1.001]);
legend('Case 1','Case 2','Case 3');

text{1} = ['CORA,ROBE22,gamma10e3,1,',num2str(compTime(1)),',',num2str(widthEnd(1))];
text{2} = ['CORA,ROBE22,gamma10e5,1,',num2str(compTime(2)),',',num2str(widthEnd(2))];
text{3} = ['CORA,ROBE22,gamma10e7,1,',num2str(compTime(3)),',',num2str(widthEnd(3))];

%------------- END OF CODE --------------