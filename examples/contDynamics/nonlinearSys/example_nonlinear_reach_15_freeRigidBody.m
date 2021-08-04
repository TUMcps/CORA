function completed = example_nonlinear_reach_15_freeRigidBody
% example_nonlinear_reach_15_freeRigidBody - example of nonlinear reachability 
% analysis; 
%
% Syntax:  
%    example_nonlinear_reach_15_freeRigidBody
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Roberto Lampariello
% Written:      ---
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

dim_x = 3;

% Model Parameters
params.tStart=0;
params.tFinal=48; %33;
% options.x0=[0; 0.1; 0];
% options.x0=[2; 4; 4; 1; 2; 0.2]*pi/180;
% params.R0=zonotope([options.x0,0.2*eye(dim)]);

% Z0{1} = zonotope([[0; 0.1; 0], 0.002*eye(dim_x)]);
% Z0{1} = zonotope([1.4 0.6 0; 2.3 0 0.05]);
params.R0 = zonotope([[0; 0.1; 0], 0.002*eye(dim_x)]);
figure; subplot(1,2,1); plot(params.R0,[1,2]); subplot(1,2,2); plot(params.R0,[1,3]);
% split initial set:
% 1st - right half in x1/x2, right upper half in x1/x3
params.R0 = zonotope([0.002; 0.1; 0.001], diag([0.000,0.002,0.001]));
params.R0 = zonotope([0.00125; 0.1; 0.001], diag([0.00075,0.002,0.001]));
subplot(1,2,1); hold on; plot(params.R0,[1,2],'r');
subplot(1,2,2); hold on; plot(params.R0,[1,3],'r');
close;
% 2nd - left half in x1/x2, left lower half in x1/x3
% params.R0 = zonotope([-0.00125; 0.1; -0.001], diag([0.00075,0.002,0.001]));
% params.R0 = polyZonotope(params.R0);

% uncertain inputs
% params.uTrans = 0;
params.U = zonotope(zeros(3,1));
dim_u = dim(params.U);

% artificial guard set version
params.startLoc = 1;


% Reachability Settings
options.timeStep = 0.025;
options.taylorTerms = 10;
options.zonotopeOrder = 200;

options.alg = 'poly';
options.tensorOrder = 3;

options.errorOrder = 5;
options.intermediateOrder = 10;

options.reductionInterval = Inf;
options.maxError = 1e-1*ones(dim_x,1); % 1e-3*ones(dim_x,1);

% settings for polynomial zonotopes
if strcmp(options.alg,'poly')
    polyZono.maxDepGenOrder = 50;
    polyZono.maxPolyZonoRatio = 0.001;
    polyZono.restructureTechnique = 'reduceFullGirard';
    options.polyZono = polyZono;
end

options.verbose = true;

% hybrid options
options.guardIntersect = 'zonoGirard';
options.enclose = {'pca'};


% Initialize System
rigid_body = nonlinearSys(@freeRigidBodyparamEq,dim_x,dim_u);


% Hybrid Automaton
startpoint = [0; 0.08; 0.01];
endpoint = [0.02; 0.12; 0.01];

c = [0.04 -0.01 0.01];
c = c / vecnorm(c',2);
d = -0.018;
% guard = halfspace(c,d);
guard = conHyperplane(c,d);
inv1 = mptPolytope(c,d);
inv2 = mptPolytope([1 0 0],1000);

reset.A = eye(dim_x); reset.b = zeros(dim_x,1);

trans{1} = transition(guard,reset,2);

loc{1} = location('loc1',inv1,trans,rigid_body);
loc{2} = location('loc2',inv2,[],rigid_body);

HA = hybridAutomaton(loc);

% Reachability Analysis
tic;
% R = reach(rigid_body, params, options);
[R,res] = reach(HA,params,options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);

% options.alg = [options.alg '-adaptive'];
% tic;
% Radap = reach(rigid_body, params, options);
% tComp = toc;
% disp(['computation time of reachable set: ',num2str(tComp)]);

% Simulations
% params.tFinal = 100;
simOpt.points = 200;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 6;
simRes = simulateRandom(rigid_body, params, simOpt);


figure; hold on; box on;
plot(simRes,[1,2,3],'k');
xlabel("x"); ylabel("y"); zlabel("z");
xlim([0,0.05]); ylim([-0.05,0.1]); zlim([0,0.15]);
plot(guard,[1,2,3],'b');
close;

% Visualization
% figure;
% for plotRun=1:dim_x
%     subplot(1,dim_x,plotRun); hold on; box on;
%     % plot reachable set
% %     plotOverTime(Radap,plotRun,'m','EdgeColor','m');
%     % plot reachable set
%     plotOverTime(R,plotRun,'b','EdgeColor','b');
%     % plot simulations
%     plotOverTime(simRes,plotRun,'y');
%     % label axes
%     xlabel("t"); ylabel("x_" + plotRun);
% end


projDimsAll = {[1,2],[2,3],[1,3]};

% figure;
% for plotRun=1:length(projDimsAll)
%     subplot(1,2,plotRun);
%     projDims = projDimsAll{plotRun};
%     plot(simRes,projDims,'k');
%     for i=1:length(R.timeInterval.set)
%         plot(R.timeInterval.set{i},projDims,'y');
%         pause(0.1);
%     end
% end

figure;
for plotRun=1:length(projDimsAll)
    subplot(1,length(projDimsAll),plotRun); hold on; box on;
    
    % projection
    projDims = projDimsAll{plotRun};
    
    % plot reachable sets
%     plot(Radap,projDims,'m','Filled',false);
    
    % plot reachable sets
    plot(R,projDims,'b','Filled',false);
    
    % plot initial set
    plot(params.R0,projDims,'w','Filled',true,'EdgeColor','k');
    
    % plot simulation results      
    plot(simRes,projDims,'k');
    
    % plot guard set
%     plot(guard,projDims,'r');

    % label axes
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end


%example completed
completed = true;

%------------- END OF CODE --------------

