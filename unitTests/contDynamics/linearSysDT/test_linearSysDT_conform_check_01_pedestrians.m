function res = test_linearSysDT_conform_check_01_pedestrians
% test_linearSysDT_conform_check_01_pedestrians - unit test function for 
%   conformance checking of linear discrete-time systems.
%
%   Checks whether recorded movement of pedestrians is within the reachable
%   set as described in [1].
%
% Syntax:
%    res = test_linearSysDT_conform_check_01_pedestrians
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] S. B. Liu, H. Roehm, C. Heinzemann, I. Lütkebohle, J. Oehlerking 
%        and M. Althoff, "Provably safe motion of mobile robots in human 
%        environments," 2017 IEEE/RSJ International Conference on 
%        Intelligent Robots and Systems (IROS), 2017, pp. 1351-1357.

% Authors:       Matthias Althoff
% Written:       07-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
path = [CORAROOT filesep 'models' filesep 'testCases' filesep 'pedestrians'];

% load test suite
load([path filesep 'Pellegrini2009Test'],'Pellegrini2009Test');

% set maximum acceleration
a_max = 2.5; % m/s^2

%% Conformance settings
options.confAlg = 'dyn';
options.zonotopeOrder = 200;
options.norm = 'interval';
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
params.W = E*W; % disturbance for 1 m/s^2 uncertain acceleration
params.V = zonotope([zeros(2,1),0.1*eye(2)]); % measurement uncertainty
params.R0conf = zonotope(zeros(4,1)); % initial state uncertainty
    
%% Conformance checking
% the result should be true according to the unit test
% res = testLong_linearSysDT_conform_02_pedestrians_bruteForce
res = conform(pedestrian,'check',params,options); 

end
        
% ------------------------------ END OF CODE ------------------------------
