function completed = example_nonlinear_conform_01_autonomousCar
% example_nonlinear_conform_01_autonomousCar - example of 
%     nonlinear conformance synthesis for following a reference trajectory; 
%     this example is also a unit test function.
%
%     This example is taken from [1].
%
% Syntax:
%    completed = example_nonlinear_conform_01_autonomousCar
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] M. Althoff and J. M. Dolan. Reachability computation of low-order 
%        models for the safety verification of high-order road vehicle 
%        models. In Proc. of the American Control Conference, 
%        page 3559–3566, 2012.

% Authors:       Matthias Althoff
% Written:       28-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Warning
disp('Warning: This example typically runs for many hours. Please comment the following lines to run it.');
completed = false;
return

% set path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'autonomousDriving'];

% load test suite
load([path filesep 'ACC2012Test'],'ACC2012Test');

% intial set
params.R0 = zonotope([[0; 0; 0; 15; 0; 0; 0], ...
    diag([0.2, 0.2, 0, 0.2, 0.05, 0.05, 0.02])]);

% uncertain inputs (disturbance W and sensor noise V are included in U 
% according to the definition of u in the model)
W = zonotope([[0; 0; 0; 0; 0; 0; 0], diag([0, 0, 0, 0, 0, 0, 0])]); % disturbance
V = zonotope([zeros(5,1), 4*diag([0.02, 0.02, 0.05*pi/180, 0.05*pi/180, 0.02])]); % sensor noise
U_tmp = cartProd(zonotope(zeros(5,1)), V);
params.U = cartProd(U_tmp, W);

%% Reachability settings

options.timeStepDivider = 2; % has to be one or greater
options.taylorTerms = 4;
options.zonotopeOrder = 800;
options.reductionTechnique = 'girard';
options.errorOrder = 2;
options.intermediateOrder = 10;
options.postProcessingOrder = 1;
%options.lagrangeRem.simplify = 'optimize';

options.reachAlg = 'lin';
options.tensorOrder = 3;

%% Conformance settings
options.confAlg = 'RRT';
options.U_increment.val = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1]; % increments for additive disturbance
options.U_increment.ind = 11:17; % maps state indices to input indices to achieve reachset conformance 
params.testSuite = ACC2012Test;


%% RRT settings
%options.points = 100;
options.points = 2;
options.vertSamp = true;
options.stretchFac = 2;
options.convertFromAbstractState = @init_MB_BMW; % handle to mFile to convert abstract states to full states (requires assumptions)
options.convertToAbstractState = @MB_to_ST; % handle to mFile to convert full states to abstract states
path = [CORAROOT filesep 'unitTests' filesep 'contDynamics' filesep 'nonlinearSys' filesep 'data'];
options.preComputedRRT = [path filesep 'simResRRT_noAdditionalInformation_27June2023'];

%% System dynamics
vehModelACC2012 = nonlinearSys(ACC2012Test{1}.model,7,17); % abstract model
CommonRoadMB2 = nonlinearSys(@vehicleDynamics_MB_controlled_BMW,29,10); % reference model
options.refModel = CommonRoadMB2; 
    
%% Conformance checking
tic
[res, R, simRes] = conform(vehModelACC2012,params,options); % RRT results should also be returned
tComp = toc

if res
    disp('Model is reachset conformant');
end


%% Visualization 
refDims = {[1 2],[3 4],[4 5]};
dims = {[1 2],[5 6],[6 4]};

%dims = {[1 2],[3 4],[5 6]};

% loop through test cases
for iCase = 1:length(params.testSuite)
    % reference trajectory
    xRef = params.testSuite{iCase}.u;
    % loop through projected dimensions
    for k = 1:length(dims)

        figure; hold on; box on
        projDim = dims{k};
        projRefDim = refDims{k};

        % plot reachable sets
        useCORAcolors("CORA:contDynamics")
        plot(R{iCase},projDim,'DisplayName','Reachable set');

        % plot simulation results      
        plot(simRes{iCase},projDim,'ko','LineStyle','none','DisplayName','Simulations');
        
        % plot initial set
        plot(R{iCase}(1).R0,projDim, 'DisplayName','Initial set');
        
        % plot reference trajectory
        plot(xRef(:,projRefDim(1)),xRef(:,projRefDim(2)),'r', ...
            'LineWidth',2,'DisplayName','Reference trajectory');
        
        % label plot
        xlabel(['x_{',num2str(projDim(1)),'}']);
        ylabel(['x_{',num2str(projDim(2)),'}']);
        legend()
    end
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
