function completed = example_linearSysDT_conform_01_pedestrians
% example_linearSysDT_conform_01_pedestrians - example of 
%     linear conformance synthesis of pedestrians; 
%     this example is also a unit test function.
%
%     This example is inspired by [1].
%
% Syntax:
%    completed = example_linearSysDT_conform_01_pedestrians
%
% Inputs:
%    -
%
% Outputs:
%    completed - boolean 
%
% References:
%    [1] S. B. Liu, H. Roehm, C. Heinzemann, I. Lütkebohle, J. Oehlerking 
%        and M. Althoff, "Provably safe motion of mobile robots in human 
%        environments," 2017 IEEE/RSJ International Conference on 
%        Intelligent Robots and Systems (IROS), 2017, pp. 1351-1357.
%    [2] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.

% Authors:       Matthias Althoff
% Written:       29-June-2023
% Last update:   21-March-2024 (LL, adapt to restructuring of "conform")
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'pedestrians'];

% load test suite
load([path filesep 'Pellegrini2009Test'],'Pellegrini2009Test');

%% Conformance settings
params.testSuite = cell(100,1);
for m=1:100
    % add nominal inputs equal to zeros to each test case since minimal input
    % dimension for linearSysDT is 1
    params.testSuite{m} = set_u(Pellegrini2009Test{m}, zeros(size(Pellegrini2009Test{m}.y,1),1));
end


%% System dynamics
% create pedestrian model (see Lecture "Formal Methods for Cyber-Physical
% Systems - Conformance Checking")
% sample time
dt = Pellegrini2009Test{1}.sampleTime;
% discrete time system matrix from continuous system matrix Ac
Ac = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
A = expm(Ac*dt);
B = zeros(4,1);
C = [1 0 0 0; 0 1 0 0];
D = zeros(2,1);
E = [dt,0,0.5*dt^2,0;0,dt,0,0.5*dt^2;0,0,dt,0;0,0,0,dt]; % input conversion from continuous time to discrete time
F = eye(2);

% instantiate linear discrete-time dynamics
pedestrian = linearSysDT(A,B,[],C,D,[],E,F,dt);
params.tFinal = 2;

% uncertain inputs (disturbance W and sensor noise V are included in U 
% according to the definition of u in the model)
Zcircle = zonotope(ellipsoid(eye(2)),'outer:norm',20); % zonotope template (uniform acceleration limit)
params.W = cartProd(interval([0;0]),Zcircle); % disturbance template
params.V = zonotope([zeros(2,1),eye(2)]); % measurement uncertainty template
params.R0 = zonotope([zeros(4,1),eye(4)]); % initial state uncertainty template for conformance

% options to weight cost function of different time steps
maxNrOfTimeSteps = ceil(params.tFinal/dt); % maximum number of timeSteps
options.cs.w = ones(maxNrOfTimeSteps+1,1);
options.cs.constraints = 'half';
    
%% Conformance synthesis
% interval norm
options.cs.cost = 'interval';
[params_interval, results] = conform(pedestrian,params,options); 
union_y_a = results.unifiedOutputs;
% Frobenius norm
options.cs.P = eye(2);
options.cs.cost = 'frob';
params_frob = conform(pedestrian,params,options); 

% %% for debugging: check conformance
% params_interval.testSuite = params.testSuite;
% res = isconform(pedestrian,params_interval,options); 
% params_frob.testSuite = params.testSuite;
% res = isconform(pedestrian,params_frob,options);

%% Compute reachable set using obtained parameters
options.zonotopeOrder = inf;
options = rmfield(options, "cs");
params_interval = rmfield(params_interval, "testSuite");
R_interval = reach(pedestrian, params_interval, options);
params_frob = rmfield(params_frob, "testSuite");
R_frob = reach(pedestrian, params_frob, options);

%% Plot test cases and reachable sets
% select projections
dims = {[1 2]};

for k = 1:length(dims)
    figure;
    
    % create separate plot for each time step
    for iStep = 1:length(R_interval.timePoint.time)
    
        subplot(3,2,iStep); hold on; box on
        projDims = dims{k};

        % plot reachable set of interval norm
        plot(R_interval.timePoint.set{iStep},projDims);
        
        % plot reachable set of Frobenius norm
        plot(R_frob.timePoint.set{iStep},projDims,'g');

        % plot unified test cases
        plot(squeeze(union_y_a{1}(iStep,projDims(1),:)),...
            squeeze(union_y_a{1}(iStep,projDims(2),:)),'kx');

        % label plot
        xlabel(['x_{',num2str(projDims(1)),'}']);
        ylabel(['x_{',num2str(projDims(2)),'}']);
        title(['Time step ',num2str(iStep)]);
        legend('interval norm','Frobenius norm')
    end
end


% example completed
completed = true;

end

% ------------------------------ END OF CODE ------------------------------
