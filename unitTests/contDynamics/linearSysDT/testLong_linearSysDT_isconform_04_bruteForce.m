function res = testLong_linearSysDT_isconform_04_bruteForce
% testLong_linearSysDT_isconform_04_bruteForce - unit test 
%   function for conformance checking of linear discrete-time systems. We
%   compare the results against the brute force method where for each test
%   case a separate reachability analysis is conducted.
%
%   Checks whether recorded movement of pedestrians is within the reachable
%   set as described in [1].
%
% Syntax:
%    testLong_linearSysDT_isconform_04_bruteForce
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
params.testSuite = Pellegrini2009Test;

%% System dynamics
% create pedestrian model (see Lecture "Formal Methods for Cyber-Physical
% Systems - Conformance Checking")
% sample time
dt = Pellegrini2009Test(1).dt;
% discrete time system matrix from continuous system matrix Ac
Ac = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
A = expm(Ac*dt);
B = [];
C = [1 0 0 0; 0 1 0 0];
D = [];
E = [dt,0,0.5*dt^2,0;0,dt,0,0.5*dt^2;0,0,dt,0;0,0,0,dt]; % input conversion from continuous time to discrete time
F = eye(2);

% instantiate linear discrete-time dynamics
pedestrian = linearSysDT(A,B,[],C,D,[],E,F,dt);
params.tFinal = 2;

% uncertain inputs (disturbance W and sensor noise V are included in U 
% according to the definition of u in the model)
Zcircle = a_max*zonotope(ellipsoid(eye(2)),'outer:norm',20); % zonotope approximating the set of possible accelerations
W = cartProd(interval([0;0]),Zcircle); % disturbance for continuous time
params.W = W; % disturbance for 1 m/s^2 uncertain acceleration
params.V = zonotope([zeros(2,1),0.1*eye(2)]); % measurement uncertainty
params.R0 = zonotope(zeros(4,1));
    
%% Conformance checking of both Approaches for a_max = 2.5 (both methods should be conformant)
% fast approach
confAlg = 'dyn';
options.zonotopeOrder = 200;
res_fast = isconform(pedestrian,params,options,confAlg); 
% brute force approach
options.reachAlg = 'standard';
options.timeStepDivider = 1;
options.zonotopeOrder = inf;
options.postProcessingOrder = inf;
options.reductionTechnique = 'girard';
confAlg = 'BF';
res_bruteForce = isconform(pedestrian,params,options,confAlg);
% results identical?
assert(res_fast == res_bruteForce);

%% Conformance checking of both Approaches for a_max = 2.4 (both methods should not be conformant)
% set maximum acceleration
a_max = 2.4; % m/s^2
Zcircle = a_max*zonotope(ellipsoid(eye(2)),'outer:norm',20); % zonotope approximating the set of possible accelrations
W = cartProd(interval([0;0]),Zcircle); % disturbance for continuous time
params.W = W; % disturbance for 1 m/s^2 uncertain acceleration
% fast approach
confAlg = 'dyn';
res_fast = isconform(pedestrian,params,options,confAlg); 
% brute force approach
confAlg = 'BF';
res_bruteForce = isconform(pedestrian,params,options,confAlg);
% results identical?
assert(res_fast == res_bruteForce);


% test completed
res = true;
        
% ------------------------------ END OF CODE ------------------------------
