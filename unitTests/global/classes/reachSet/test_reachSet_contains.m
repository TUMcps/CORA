function res = test_reachSet_contains
% test_reachSet_contains - unit test function for contains
%
% Syntax:
%    res = test_reachSet_contains()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       28-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple continuous-time linear system
A = [0.1 1; -1 0.1];
B = 1;
sys_ct = linearSys('linearSys',A,B);

% model parameters and reachability options
params.tFinal = 1;
params.R0 = zonotope([10;5],0.2*eye(2));
options.linAlg = 'adaptive';
options.error = 0.1;

% compute reachable set
R = reach(sys_ct,params,options);

% compute simulation
simOpt.points = 25;
simRes = simulateRandom(sys_ct,params,simOpt);

% visualization
% figure; hold on; box on;
% useCORAcolors("CORA:contDynamics");
% plot(R);
% plot(R.R0);
% plot(simRes);
% close;

% check containment
tol = 1e-10;
res = contains(R,simRes,'exact',tol);


% simple discrete-time linear system
A = [0.95 -0.15;...
    0.15 0.95];
dim_x = length(A);
B = 1;
dt = 0.04;
sys_dt = linearSysDT(A,B,dt);

% model parameters and reachability options
params.tFinal = 1;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = zonotope(zeros(dim_x,1),0.02*diag([0.1, 0.3]));
clear options
options.zonotopeOrder = 10;

% compute reachable set
R = reach(sys_dt,params,options);

% compute simulation
simRes = simulateRandom(sys_dt,params,simOpt);

% visualization
% figure; hold on; box on;
% useCORAcolors("CORA:contDynamics");
% plot(R);
% plot(R.R0);
% plot(simRes);
% close;

% check containment
tol = 1e-10;
res = res & contains(R,simRes,'exact',tol);

% ------------------------------ END OF CODE ------------------------------
