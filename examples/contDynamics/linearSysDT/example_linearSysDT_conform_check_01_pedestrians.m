function completed = example_linearSysDT_conform_check_01_pedestrians
% example_linearSysDT_conform_check_01_pedestrians - example of 
%     linear conformance checking of pedestrians; 
%     this example is also a unit test function.
%
%     This example is taken from [1].
%
% Syntax:
%    completed = example_linearSysDT_conform_check_01_pedestrians
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

% Authors:       Matthias Althoff
% Written:       29-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'pedestrians'];

% load test suite
load([path filesep 'Pellegrini2009Test'],'Pellegrini2009Test');

% set maximum velocity
a_max = 2.5; % m/s^2

%% Conformance settings
options.confAlg = 'dyn';
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
Zcircle = a_max*zonotope(ellipsoid(eye(2)),20); % zonotope approximating the set of possible accelrations
W = cartProd(interval([0;0]),Zcircle); % disturbance for continuous time
params.W = E*W; % disturbance for a_max uncertain acceleration
params.V = zonotope([zeros(2,1),0.1*eye(2)]); % measurement uncertainty
params.R0conf = zonotope(zeros(4,1)); % initial state uncertainty template for conformance 

%% Conformance checking
res = conform(pedestrian,'check',params,options); 

if res
    disp('Model is reachset conformant');
end

% example completed
completed = true;
end


% ------------------------------ END OF CODE ------------------------------
