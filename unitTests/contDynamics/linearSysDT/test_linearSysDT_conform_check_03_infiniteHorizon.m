function res = test_linearSysDT_conform_check_03_infiniteHorizon
% test_linearSysDT_conform_check_03_infiniteHorizon - unit test function 
%   for conformance synthesis of linear discrete-time systems according to 
%   [1] and [2].
%
%   Conformance for an infinite horizon according to Corollary 1 in [2] is
%   checked
%
% Syntax:
%    res = test_linearSysDT_conform_check_03_infiniteHorizon
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
params.tFinal = dt*100;

% the system has known uncertainties
W = zonotope(rand(4,5)); % disturbance 
V = zonotope(rand(2,3)); % measurement uncertainty 
R0 = zonotope(rand(4,1)); % initial state uncertainty; initial state uncertainty has no effect for linear systems

% options to weight cost function of different time steps
maxNrOfTimeSteps = ceil(params.tFinal/dt); % maximum number of timeSteps

%% Create test cases
% because the system has no imaginary eigenvalues, the extreme cases
% suffice to reproduce the actual uncertainties
% combine zonotopes for R0, W, and V in a single zonotope
S = cartProd(W, V);
% compute vertices of combined zonotope
Vmat = vertices(S);
% set initial state
params.R0 = R0;
% init test cases
nrOfTests = size(Vmat,2);
synthesizedTests = cell(nrOfTests,1);
% loop over each vertex
for i = 1:nrOfTests
    % overwrite parameters
    params.W = zonotope(Vmat(1:4,i)); % disturbance 
    params.V = zonotope(Vmat(5:6,i)); % measurement uncertainty  
    % simulate system
    simOpt.points = 1;
    simRes = simulateRandom(pedestrian,params,simOpt);
    % construct input vector
    uVec = zeros(maxNrOfTimeSteps+1,1);
    % save in test case
    synthesizedTests{i} = testCase(simRes.y{1}, uVec, simRes.x{1}, dt);
end

    
%% Conformance check (result should be true)
options.confAlg = 'dyn';
options.zonotopeOrder = 200;
params.testSuite = synthesizedTests;
% store uncertainty
params.W = enlarge(W,1.001); % disturbance enlarged by 0.1%
params.V = enlarge(V,1.001); % measurement uncertainty enlarged by 0.1%
params.R0conf = R0; % initial state uncertainty 
% set time horizon to inf;
params.tFinal = inf;
% overall result
res = conform(pedestrian,'check',params,options); 

end
        
% ------------------------------ END OF CODE ------------------------------
