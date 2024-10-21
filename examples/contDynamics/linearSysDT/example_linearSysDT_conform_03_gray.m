function completed = example_linearSysDT_conform_03_gray
% example_linearSysDT_conform_03_gray - example script for reachset
%   conformance synthesis using a gray-box model. This method identifies
%   parameters of a gray-box model using fmincon in an outer loop as in 
%   Sec. III.B of [1] and updates the required uncertain sets according 
%   to [2],[3] in an inner loop. There exists a corresponding unit test.
%
% Syntax:
%    completed = example_linearSysDT_conform_03_gray
%
% Inputs:
%    None
%
% Outputs:
%    completed - boolean
%
% References:
%    [1] S. B. Liu and M. Althoff, "Online Verification of 
%        Impact-Force-Limiting Control for Physical Human-Robot 
%        Interaction," 2021 IEEE/RSJ International Conference on 
%        Intelligent Robots and Systems (IROS), 2021, pp. 777-783.
%    [2] M. Althoff. Checking and Establishing Reachset Conformance in
%        CORA 2023. In Proc. of the 10th International Workshop on 
%        Applied Verification of Continuous and Hybrid Systems, 2023.
%    [3] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.
%

% Authors:       Matthias Althoff
% Written:       30-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% System
% set seed for determinism 
rng(2023);
n_m = 20; 
n_s = 10;
n_k = 10; % number of prediction time steps for one test case

% define state space system 
dt = 0.1;

% system matrices
A = [0.5, -0.1;
      0.03, 0.4];
B = [ -0.01, -0.03;
      -0.02, 0.03];
C = [ 24, -110;
     150,  50];
D = [ 0, -2;
      2,  0];
  
% instantiate system  
sys = linearSysDT(linearSys(A,B,[],C,D), dt);

% parameters of the system
params.tStart       = 0;
params.tFinal       = dt * (n_k-1);
params.R0           = zonotope([ones(2,1),0.1*eye(2)]);
params.U            = zonotope(zeros(2, 1), eye(2));

% Initial estimations of system matrices, e.g., obtained from a system
% identification
A_0     = A + [0.1, -0.1; 0, -0.2];
B_0     = B + [0.2, -0.2; 0.1, -0.1];
C_0     = C + [0, 0; 0.1, -0.1];
D_0     = D + [-0.1, -0.1; -0.1, -0.1];
% parameter vector is simply the concatenation of entries of all system
% matrices
p_0     = [ reshape(A_0, [], 1); ...
            reshape(B_0, [], 1); ...
            reshape(C_0, [], 1); ...
            reshape(D_0, [], 1)];

% options for the system identification
dim_x = 2;
dim_u = 2;
dim_y = 2;
n_params = dim_x^2 + dim_x*dim_u + dim_y*dim_x + dim_y*dim_u;

% reachability and conformance options
options = struct();
options.zonotopeOrder = 50; %zonotope order
options.tensorOrder = 3;
options.tensorOrderOutput = 3;
options.errorOrder = 3;
options.cs.cost = 'interval';
options.cs.constraints = 'half';
options.cs.timeout = 300; %timeout for fmincon when estimating c
options.cs.set_p = @(p,params) aux_set_p(p, params, [2 2 2], dt);
options.cs.p0 = p_0;
options.p_min               = -100*ones(n_params, 1);
options.p_max               = 100*ones(n_params, 1);
options.p_min(1:dim_x^2)    = -1*ones(dim_x^2 , 1); % bounds for system matrices are chosen tighter
options.p_max(1:dim_x^2)    = 1*ones(dim_x^2 , 1); % bounds for system matrices are chosen tighter
confAlg = 'graySeq';

% options to weight cost function of different time steps
options.cs.w = ones(1,n_k);

%% Create test cases
% options for creating the test data
options_testS.p_extr = 0.4;
options_testS.p_biased = 0;
options_testS.stateSet = interval(-10 * ones(2,1), 10 * ones(2,1));
params.testSuite = createTestSuite(sys, params, n_k, n_m, ...
        n_s, options_testS);

%% Conformance synthesis
% Perform conformance synthesis for gray-box model
[params_gray, results] = conform(sys,params,options,confAlg);
sys_gray = results.sys; 

% compute output deviation y_a
testCase_ya = params.testSuite{1}.compute_ya(sys_gray);
y_a = testCase_ya.y_a;

%% Compute reachable set using obtained parameters
options.zonotopeOrder = inf;

% compute reachable set of output deviations y_a
R_gray = reach(sys_gray, params_gray, options);

%% Plot results
% select projections
dims = {[1 2]};

for k = 1:length(dims)
    % create separate plot for each time step
    figure;

    for iStep = 1:min(length(R_gray.timePoint.time),6)
    
        subplot(3,2,iStep); hold on; box on
        projDims = dims{k};

        % plot reachable set
        plot(R_gray.timePoint.set{iStep},projDims, 'Color','r');
        
        % plot unified outputs
        plot(squeeze(y_a(iStep,projDims(1),:)),...
            squeeze(y_a(iStep,projDims(2),:)),'Marker','.','LineStyle', 'none');

        % label plot
        xlabel(['x_{',num2str(projDims(1)),'}']);
        ylabel(['x_{',num2str(projDims(2)),'}']);
        title(['Time step ',num2str(iStep)]);
    end
end

% Example completed
completed = true;

end


% Auxiliary functions -----------------------------------------------------

function [sys,params] = aux_set_p(p, params, dims, dt)
[A, B, C, D] = getMatricesFromP_twoDimExample(p,dims);
sys = linearSysDT(A,B,[],C,D,[],dt);
end

function synthesizedTests = aux_createTestCases(sys, params)
params = rmfield(params, 'getMatricesFromP');
    % because the system has no imaginary eigenvalues, the extreme cases
    % suffice to reproduce the actual uncertainties
    % combine zonotopes for R0, W, and V in a single zonotope
    
    % obtain dimensions of x, w, v
    dim_x = dim(params.R0);
    dim_u = size(sys.B,2);
    dim_w = dim(params.W);
    dim_v = dim(params.V);
    % Cartesian product of initial set, disturbance set, and measurement
    % uncertainty
    S_tmp = cartProd(params.R0, params.W);
    S = cartProd(S_tmp, params.V);
    % compute vertices of combined zonotope
    Vmat = vertices(S);
    % nr of tests 
    nrOfTests = size(Vmat,2);
    % init result
    synthesizedTests = cell(nrOfTests,1);
    % maximum number of timeSteps
    maxNrOfTimeSteps = ceil(params.tFinal/sys.dt); 
    % loop over each vertex
    for i = 1:nrOfTests
        % overwrite parameters
        params.R0 = zonotope(Vmat(1:dim_x,i)); % initial state
        params.W = zonotope(Vmat(dim_x+1:dim_x+dim_w,i)); % disturbance 
        params.V = zonotope(Vmat(dim_x+dim_w+1:dim_x+dim_w+dim_v,i)); % measurement uncertainty  
        % simulate system
        simOpt.points = 1;
        simRes = simulateRandom(sys, params, simOpt);
        % construct input vector
        uVec = zeros(maxNrOfTimeSteps+1,dim_u);
        % save in test case
        synthesizedTests{i} = testCase(simRes.y{1}, uVec, simRes.x{1}, sys.dt);
    end
end

% ------------------------------ END OF CODE ------------------------------
