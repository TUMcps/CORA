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
% Last update:   25-March-2024 (LL, adapt to new conform function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% initialize partial results
resPartial = [];
rng(1)

%% System dynamics
% create pedestrian model (see Lecture "Formal Methods for Cyber-Physical
% Systems - Conformance Checking")
% sample time
dt = 0.4;
% discrete time system matrix from continuous system matrix Ac
Ac = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
A = expm(Ac*dt);
B = [eye(4) zeros(4,2)];
C = [1 0 0 0; 0 1 0 0];
D = [zeros(2,4) eye(2)];

% instantiate linear discrete-time dynamics
pedestrian = linearSysDT(A,B,[],C,D,dt);
params.tFinal = 2;

% the system has no uncertainty
W = zonotope(10*randn(4,1)); % process disturbance 
V = zonotope(10*rand(2,1)); % measurement uncertainty 
U = cartProd(W,V); % combined ditrurbances
R0 = zonotope(10*randn(4,1)); % initial state uncertainty 

% store fixed disturbance, measurement error, and initial state in params
params.U = U; % disturbance
params.R0 = R0; % initial state uncertainty 
params.tStart = 0;

% options to weight cost function of different time steps
maxNrOfTimeSteps = ceil(params.tFinal/dt); % maximum number of timeSteps
options.cs.w = ones(maxNrOfTimeSteps+1,1);
options.cs.constraints = 'half';

%% Create test case
% simulate system
simOpt.points = 1;
traj = simulateRandom(pedestrian,params,simOpt);
% construct input vector
uVec = zeros(maxNrOfTimeSteps+1,6);
traj = trajectory(uVec, traj.x, traj.y, [], dt);

%% Conformance synthesis (default: interval norm)
options.zonotopeOrder = 200;
params.testSuite = traj;

% add uncertainty templates
W_guess = zonotope(rand(4,10)); % disturbance 
V_guess = zonotope(rand(2,10)); % measurement uncertainty 
params.U = cartProd(W_guess, V_guess);
params.R0 = zonotope([rand(4,10)]); % initial state uncertainty 
tic
params_new = conform(pedestrian,params,options);

%% check whether all parameter sets have no uncertainty up to some floating point errors
accuracy = 1e-8;

% check whether there is uncertainty in the identified R0 and U
assert((max(rad(interval(params_new.R0))) < accuracy));% R0
assert((max(rad(interval(params_new.U))) < accuracy));% U

% check whether centers match up to accuracy (must not be true)
%resPartial(end+1) = (max(abs(center(params_new.R0))) < accuracy); % R0 is relative to x_0
%resPartial(end+1) = (max(abs(center(params_new.U) - center(U))) < accuracy); % U

% check whether the measurements can be reconstructed from the
% identification restults
params_new.R0 = params_new.R0 + traj.x(:,1,1); 
params_new = rmfield(params_new, 'testSuite');
traj_new = simulateRandom(pedestrian,params_new,simOpt);
assert((max(abs(traj_new.y - traj.y),[],'all') < accuracy));

% overall result
res = all(resPartial);

end
        
% ------------------------------ END OF CODE ------------------------------
