function res = test_nonlinearSysDT_reach_01_cstrDisc
% test_nonlinearSysDT_reach_cstrDisc - unit test for nonlinear discrete time
%    reachability analysis from [1, Sec.6]. 
%
% Syntax:  
%    res = test_nonlinearSysDT_reach_cstr
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Reference:
%    [1] J.M. Bravo, Robust MPC of constrained discrete-time
%        nonlinear systems based on approximated reachable sets, 2006.

% Author:       Niklas Kochdumper, Matthias Althoff
% Written:      30-January-2018
% Last update:  20-March-2020 (MA, simulateRandomDT from inherited class)
%               23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameters --------------------------------------------------------------

params.tFinal = 0.03;
params.R0 = zonotope([-0.15;-45],diag([0.005;3]));
params.U = zonotope(zeros(2,1),diag([0.1;2]));


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 100;
options.errorOrder = 5;
options.tensorOrder = 3;


% System Dynamics ---------------------------------------------------------

% sampling time
dt = 0.015;

fun = @(x,u) cstrDiscr(x,u,dt);

sysDisc = nonlinearSysDT('stirredTankReactor',fun,0.015);


% Reachability Analysis ---------------------------------------------------

R = reach(sysDisc,params,options);


% Simulation --------------------------------------------------------------

simOpt.points = 50;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 0;

simRes = simulateRandom(sysDisc, params, simOpt);


% Verification ------------------------------------------------------------

res = true;

for i = 1:simOpt.points
   for j = 1:length(R.timePoint.set)
      if ~in(R.timePoint.set{j},simRes.x{i}(j,:)')
          res = false;
          path = pathFailedTests(mfilename());
          save(path,'simRes');       
          error('Test failed!');
      end
   end
end

%------------- END OF CODE --------------