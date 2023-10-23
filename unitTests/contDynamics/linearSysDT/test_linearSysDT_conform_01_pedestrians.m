function res = test_linearSysDT_conform_01_pedestrians
% test_linearSysDT_conform_01_pedestrians - unit test function for 
%   conformance synthesis of linear discrete-time systems according to [1].
%
%   Checks whether synthesized initial set, disturbance set, and 
%   measurement noise set result in reachable sets that contain all 
%   measurements. It isfurther checked whether a counterexample is found 
%   when the synthesized sets are reduced.
%
% Syntax:
%    res = test_linearSysDT_conform_01_pedestrians
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
% Written:       11-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize partial results
resPartial = [];

% set path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'pedestrians'];

% load test suite
load([path filesep 'Pellegrini2009Test'],'Pellegrini2009Test');

%% Conformance settings
options.confAlg = 'dyn';
options.norm = 'interval';
options.zonotopeOrder = 200;
params.testSuite = Pellegrini2009Test;


%% System dynamics
% create pedestrian model (see Lecture "Formal Methods for Cyber-Physical
% Systems - Conformance Checking")
% sample time
dt = Pellegrini2009Test{1}.sampleTime;
% discrete time system matrix from continuous system matrix Ac
Ac = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
A = expm(Ac*dt);
B = [];
C = [1 0 0 0; 0 1 0 0];
D = [];

% instantiate linear discrete-time dynamics
pedestrian = linearSysDT(A,B,[],C,D,dt);
params.tFinal = 2;

% uncertain inputs (disturbance W and sensor noise V are included in U 
% according to the definition of u in the model)
E = [dt,0,0.5*dt^2,0;0,dt,0,0.5*dt^2;0,0,dt,0;0,0,0,dt]; % input conversion from continuous time to discrete time
Zcircle = zonotope(ellipsoid(eye(2)),20); % zonotope template (uniform acceleration limit)
W = cartProd(interval([0;0]),Zcircle); % disturbance template
params.W = E*W; % disturbance template for discrete time
params.V = zonotope([zeros(2,1),eye(2)]); % measurement uncertainty template
params.R0conf = zonotope([zeros(4,1),eye(4)]); % initial state uncertainty template

% options to weight cost function of different time steps
maxNrOfTimeSteps = ceil(params.tFinal/dt); % maximum number of timeSteps
options.w = ones(maxNrOfTimeSteps+1,1);
    
%% Conformance synthesis
params_interval = conform(pedestrian,params,options); 

%% Perform conformance check on obtained parameters (result should be true)
options.confAlg = 'dyn';
Vorig = params_interval.V;
params_interval.V = enlarge(Vorig,1.001); % slightly enlarge V for numerical robustness
resPartial(end+1) = conform(pedestrian,'check',params_interval,options);

%% Reduce measurement uncertainty to 99.9% (result should be false)
params_interval.V = enlarge(Vorig,0.999);
if isequal(params_interval.V, Vorig)
    params_interval.W = enlarge(params_interval.W,0.999);
end
resPartial(end+1) = ~conform(pedestrian,'check',params_interval,options);

%% Conformance synthesis using Frobenius norm
options.confAlg = 'dyn';
options.norm = 'frob';
options.P = eye(2);
params_frob = conform(pedestrian,params,options); 

%% Perform conformance check on obtained parameters (result should be true)
options.confAlg = 'dyn';
Vorig = params_frob.V;
params_frob.V = enlarge(Vorig,1.001); % slightly enlarge V for numerical robustness
resPartial(end+1) = conform(pedestrian,'check',params_frob,options);

%% Reduce measurement uncertainty to 99.9% (result should be false)
params_frob.V = enlarge(Vorig,0.999);
if isequal(params_frob.V, Vorig)
    params_frob.W = enlarge(params_frob.W,0.999);
end
resPartial(end+1) = ~conform(pedestrian,'check',params_frob,options);

% overall result
res = all(resPartial);

end
        
% ------------------------------ END OF CODE ------------------------------
