function res = test_linearSysDT_conform_02_noUncertainty
% test_linearSysDT_conform_02_noUncertainty - unit test function 
%   for conformance synthesis of linear discrete-time systems according to 
%   [1].
%
%   This is a sanity check. A single test case without uncertainty is
%   produced. The conformnace synthesis should only produce uncertain sets
%   due to floating point errors. Interestingly this is only possible for 
%   the interval norm and not the Frobenius norm. The Frobenius norm 
%   produces more conservative results when there is no uncetainty.
%
% Syntax:
%    res = test_linearSysDT_conform_02_noUncertainty
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

% the system has no uncertainty
W = zonotope(rand(4,1)-0.5*ones(4,1)); % disturbance 
V = zonotope(rand(2,1)-0.5*ones(2,1)); % measurement uncertainty 
R0 = zonotope(rand(4,1)-0.5*ones(4,1)); % initial state uncertainty 

% store fixed disturbance, measurement error, and initial state in params
params.W = W; % disturbance 
params.V = V; % measurement uncertainty 
params.R0 = R0; % initial state uncertainty 

% options to weight cost function of different time steps
maxNrOfTimeSteps = ceil(params.tFinal/dt); % maximum number of timeSteps
options.w = ones(maxNrOfTimeSteps+1,1);

%% Create test case
% simulate system
simOpt.points = 1;
simRes = simulateRandom(pedestrian,params,simOpt);
% construct input vector
uVec = zeros(maxNrOfTimeSteps+1,1);
% save in test case
sanityCheck{1} = testCase(simRes.y{1}, uVec, simRes.x{1}, dt);
    
%% Conformance synthesis (default: interval norm)
options.confAlg = 'dyn';
options.zonotopeOrder = 200;
params.testSuite = sanityCheck;
% add options for uncertainty templates
params.W = zonotope(rand(4,10)); % disturbance 
params.V = zonotope(rand(2,10)); % measurement uncertainty 
params.R0conf = zonotope([rand(4,10)]); % initial state uncertainty 
params_interval = conform(pedestrian,params,options); 

%% check whether all parameter sets have no uncertainty up to some floating point errors
accuracy = 1e-8;
resPartial(end+1) = (max(rad(interval(params_interval.R0conf))) < accuracy); % R0
resPartial(end+1) = (max(rad(interval(params_interval.W))) < accuracy); % W
resPartial(end+1) = (max(rad(interval(params_interval.V))) < accuracy); % V

% check whether centers match up to accuracy
resPartial(end+1) = (max(abs(center(params_interval.R0conf))) < accuracy); % R0conf is relative to x_0
resPartial(end+1) = (max(abs(center(params_interval.W) - center(W))) < accuracy); % W
resPartial(end+1) = (max(abs(center(params_interval.V) - center(V))) < accuracy); % V

% overall result
res = all(resPartial);

end
        
% ------------------------------ END OF CODE ------------------------------
