function text = benchmark_nonlinear_reach_ARCH23_Robertson()
% benchmark_nonlinear_reach_ARCH23_Robertson - Robertson chemical reaction
%    example from the ARCH23 competition
%
% Syntax:
%    benchmark_nonlinear_reach_ARCH23_Robertson
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%

% Authors:       Mark Wetzlinger
% Written:       17-May-2021
% Last update:   25-March-2023 (use adaptive algorithm for all cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

tFinal = 40;
params.tFinal = tFinal; 
params.U = zonotope(0);
params.R0 = zonotope([1;0;0]);

dim_x = dim(params.R0);
dim_u = dim(params.U);

options.alg = 'lin-adaptive';
options.verbose = true;

% System Dynamics ---------------------------------------------------------

sys{1} = nonlinearSys('Robertson_case1_ARCH23',@Robertson_case1,dim_x,dim_u);
sys{2} = nonlinearSys('Robertson_case2_ARCH23',@Robertson_case2,dim_x,dim_u);
sys{3} = nonlinearSys('Robertson_case3_ARCH23',@Robertson_case3,dim_x,dim_u);

nrSys = length(sys);


% Reachability Analysis ---------------------------------------------------

R = cell(nrSys,1);
compTime = zeros(nrSys,1);
widthSum = cell(nrSys,1);
widthEnd = zeros(nrSys,1);
noSets = zeros(nrSys,1);

for i=1:nrSys
    
    tic;
    R{i} = reach(sys{i}, params, options);
    compTime(i) = toc;
    
    % compute sum of concentrations of all time-interval sets
	noSets(i) = length(R{i}.timeInterval.set);

    widthSum{i} = cell(noSets(i),1);
    for iSet = 1:noSets(i)
        % width (interval) of each set
        widthSum{i}{iSet,1} = sum(interval(R{i}.timeInterval.set{iSet}));
    end
    widthEnd(i,1) = 2*rad(widthSum{i}{end});
    
end


% Visualization -----------------------------------------------------------

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
xlabel('$t$');
ylabel('$s$');
% axis([0,params.tFinal,0.999,1.001]);
legend('Case 1','Case 2','Case 3');

text{1} = ['ROBE21,gamma10e3,1,',num2str(compTime(1)),',',num2str(widthEnd(1)),',',noSets(1)];
text{2} = ['ROBE21,gamma10e5,1,',num2str(compTime(2)),',',num2str(widthEnd(2)),',',noSets(2)];
text{3} = ['ROBE21,gamma10e7,1,',num2str(compTime(3)),',',num2str(widthEnd(3)),',',noSets(3)];

% ------------------------------ END OF CODE ------------------------------
