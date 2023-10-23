function res = test_linearSysDT_conform_03_knownUncertainty
% test_linearSysDT_conform_03_knownUncertainty - unit test function 
%   for conformance synthesis of linear discrete-time systems according to 
%   [1].
%
%   This is a sanity check. Test cases with known uncertainty are
%   produced. The conformnace synthesis should reproduce the known 
%   uncertainty. Interestingly this is only possible for the interval norm 
%   and not the Frobenius norm. The Frobenius norm produces more 
%   conservative results so that the known uncertainty is no exactly 
%   reproduced.
%
% Syntax:
%    res = test_linearSysDT_conform_03_knownUncertainty
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Althoff. Checking and Establishing Reachset Conformance in
%        CORA 2023. In Proc. of the 10th International Workshop on 
%        Applied Verification of Continuous and Hybrid Systems, 2023.
%    [2] Liu et al., "Guarantees for Real Robotic Systems: Unifying Formal
%        Controller Synthesis and Reachset-Conformant Identification", 202x.

% Authors:       Matthias Althoff
% Written:       12-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize partial results
resPartial = [];

%% System dynamics
% create pedestrian model (see Lecture "Formal Methods for Cyber-Physical
% Systems - Conformance Checking")
% sample time
dt = 0.4;
% discrete time system matrix from continuous system matrix Ac
Ac = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
A = expm(Ac*dt);
B = [];
C = [1 0 0 0; 0 1 0 0];
D = [];

% instantiate linear discrete-time dynamics
pedestrian = linearSysDT(A,B,[],C,D,dt);
params.tFinal = 2;

% options for conformance synthesis
options.confAlg = 'dyn';
maxNrOfTimeSteps = ceil(params.tFinal/dt); % maximum number of timeSteps
options.w = ones(maxNrOfTimeSteps+1,1);

%% CASE 1: only disturbance
% the system has known uncertainties
W = zonotope(rand(4,5)); % disturbance 
V = zonotope(rand(2,1)); % measurement uncertainty 
R0 = zonotope(rand(4,1)); % initial state uncertainty; initial state uncertainty has no effect for linear systems

% Create test cases:
% because the system has no imaginary eigenvalues, the extreme cases
% suffice to reproduce the actual uncertainties
synthesizedTests = aux_createTestCases(pedestrian, R0, W, V, params, dt);

%% Conformance synthesis
% load params struct
params.testSuite = synthesizedTests;
params.W = W; % disturbance 
params.V = V + zonotope([zeros(2,1),eye(2)]); % measurement uncertainty 
params.R0conf = R0; % initial state uncertainty for conformance

%% Perform conformance synthesis using interval norm
params_interval = conform(pedestrian,params,options); 
% check whether synthesized parameters enclose true uncertainties
resPartial = aux_checkCorrectness(resPartial,params_interval,W,V);

%% CASE 2: only measurement uncertainty
% the system has known uncertainties
W = zonotope(rand(4,1)); % disturbance 
V = zonotope(rand(2,3)); % measurement uncertainty 
R0 = zonotope(rand(4,1)); % initial state uncertainty; initial state uncertainty has no effect for linear systems

% Create test cases:
% because the system has no imaginary eigenvalues, the extreme cases
% suffice to reproduce the actual uncertainties
synthesizedTests = aux_createTestCases(pedestrian, R0, W, V, params, dt);

%% Conformance synthesis
% load params struct
params.testSuite = synthesizedTests;
params.W = W; % + zonotope([zeros(4,1),eye(4)]); % disturbance 
params.V = V; % measurement uncertainty 
params.R0conf = R0; % initial state uncertainty for conformance
% perform conformance synthesis
params = conform(pedestrian,params,options); 

% check whether synthesized parameters enclose true uncertainties
resPartial = aux_checkCorrectness(resPartial,params,W,V);


%% CASE 3: disturbance and measurement uncertainty
% the system has known uncertainties
W = zonotope(rand(4,5)); % disturbance 
V = zonotope(rand(2,3)); % measurement uncertainty 
R0 = zonotope(rand(4,1)); % initial state uncertainty; initial state uncertainty has no effect for linear systems

% Create test cases:
% because the system has no imaginary eigenvalues, the extreme cases
% suffice to reproduce the actual uncertainties
synthesizedTests = aux_createTestCases(pedestrian, R0, W, V, params, dt);

%% Conformance synthesis
% load params struct
params.testSuite = synthesizedTests;
params.W = W; % disturbance 
params.V = V; % measurement uncertainty; add small unecrtainty for numerical stability
params.R0conf = R0; % initial state uncertainty for conformance
% perform conformance synthesis
params = conform(pedestrian,params,options); 

% check whether synthesized parameters enclose true uncertainties
resPartial = aux_checkCorrectness(resPartial,params,W,V);


%% CASE 4: initial state, disturbance, and measurement uncertainty
% the system has known uncertainties
W = zonotope(rand(4,3)); % disturbance 
V = zonotope(rand(2,3)); % measurement uncertainty 
R0 = zonotope(rand(4,3)); % initial state uncertainty; initial state uncertainty has no effect for linear systems

% Create test cases:
% because the system has no imaginary eigenvalues, the extreme cases
% suffice to reproduce the actual uncertainties
synthesizedTests = aux_createTestCases(pedestrian, R0, W, V, params, dt);

%% Conformance synthesis
% load params struct
params.testSuite = synthesizedTests;
params.W = W; % disturbance 
params.V = V; % measurement uncertainty; add small unecrtainty for numerical stability
params.R0conf = R0; % initial state uncertainty for conformance
% perform conformance synthesis
params = conform(pedestrian,params,options); 

% check whether synthesized parameters enclose true uncertainties
resPartial = aux_checkCorrectness(resPartial,params,W,V);


% overall result
res = all(resPartial);

end


% Auxiliary functions -----------------------------------------------------

function synthesizedTests = aux_createTestCases(sys, R0, W, V, params, dt)
    % because the system has no imaginary eigenvalues, the extreme cases
    % suffice to reproduce the actual uncertainties
    % combine zonotopes for R0, W, and V in a single zonotope
    S_tmp = cartProd(R0, W);
    S = cartProd(S_tmp, V);
    % compute vertices of combined zonotope
    Vmat = vertices(S);
    % nr of tests 
    nrOfTests = size(Vmat,2);
    % init result
    synthesizedTests = cell(nrOfTests,1);
    % maximum number of timeSteps
    maxNrOfTimeSteps = ceil(params.tFinal/dt); 
    % loop over each vertex
    for i = 1:nrOfTests
        % overwrite parameters
        params.R0 = zonotope(Vmat(1:4,i)); % initial state
        params.W = zonotope(Vmat(5:8,i)); % disturbance 
        params.V = zonotope(Vmat(9:10,i)); % measurement uncertainty  
        % simulate system
        simOpt.points = 1;
        simRes = simulateRandom(sys, params, simOpt);
        % construct input vector
        uVec = zeros(maxNrOfTimeSteps+1,1);
        % save in test case
        synthesizedTests{i} = testCase(simRes.y{1}, uVec, simRes.x{1}, dt);
    end
end

function resPartial = aux_checkCorrectness(resPartial,params,W,V)
    % check whether synthesized parameters enclose true uncertainties
    accuracy = 1e-6;
    resPartial(end+1) = contains(params.R0conf + zonotope([zeros(4,1),accuracy*eye(4)]),zeros(dim(params.R0conf),1)); % R0conf should contain the origin so that R0+R0condf contains R0
    resPartial(end+1) = contains(params.W + zonotope([zeros(4,1),accuracy*eye(4)]),W); % W
    resPartial(end+1) = contains(params.V + zonotope([zeros(2,1),accuracy*eye(2)]),V); % V
end
        
% ------------------------------ END OF CODE ------------------------------
