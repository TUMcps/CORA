function text = benchmark_hybrid_reach_ARCH23_brake_BRKNC01()
% benchmark_hybrid_reach_ARCH23_brake_BRKNC01 -  brake benchmark from the
%    2023 ARCH competition (non-deterministic switching)
%
% Syntax:
%    benchmark_hybrid_reach_ARCH23_brake_BRKNC01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%

% Authors:       Niklas Kochdumper
% Written:       27-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Dynamic System ----------------------------------------------------------

% parameter values
L = 0.001;
K_P = 10000;
K_I = 1000;
R = 0.5;
K = 0.02;
d_rot = 0.1;
i = 113.1167;
T_sample = 0.0001;
x0 = 0.05;
t_min = -1e-8;
t_max = 1e-7;

% system matrices
A = [-(R+K^2/d_rot)/L, 0, K_P/L, K_I/L, 0; ....
     K/(i*d_rot), 0, 0, 0, 0; ...
     0, 0, 0, 0, 0; ...
     0, 0, 0, 0, 0; ...
     0, 0, 0, 0, 0];
 
B = [0;0;0;0;0];
c = [0;0;0;0;1];

% reset
reset.A = [1 0 0 0 0; 0 1 0 0 0; 0 -1 0 0 0;0 -T_sample 0 1 0;0 0 0 0 1];
reset.c = [0; 0; x0; T_sample*x0; -T_sample];

% invariant set
inv = polytope([0 0 0 0 1],T_sample + t_max);

% guard set
guard = polytope([0 0 0 0 1; 0 0 0 0 -1], ...
                    [T_sample+t_max;-(T_sample+t_min)]);

% linear system object
sys = linearSys(A,B,c);

% hybrid automaton
trans(1) = transition(guard,reset,1);
loc(1) = location(inv,trans,sys);

HA = hybridAutomaton(loc);


% Parameter ---------------------------------------------------------------

params.tFinal = 0.1;
params.R0 = zonotope([0;0;0;0;0],diag(1e-8*[1 1 1 1 1]));
params.startLoc = 1;


% Settings ----------------------------------------------------------------

options.timeStep = T_sample/5;
options.taylorTerms = 5;
options.zonotopeOrder = 20;

options.guardIntersect = 'conZonotope';
options.enclose = {'box','pca'};
options.guardOrder = 5;


% Specification -----------------------------------------------------------

hs = halfspace([0 -1 0 0 0],-x0);

spec = specification(hs,'unsafeSet');


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = reach(HA,params,options);
tComp = toc(clock);


% Verification ------------------------------------------------------------

clock = tic;
res = true;

list = query(R,'reachSet');

for i = 1:length(list)
   temp = supremum(interval(project(list{i},2)));
   if temp > x0
      res = false;
      break;
   end
end
tVer = toc(clock);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tVer + tComp)]);
 

% Visualization -----------------------------------------------------------

figure; hold on; box on;
useCORAcolors('CORA:contDynamics')

% plot reachable set
plotOverTime(R,2);

% plot specification
spec = specification(interval([0; x0], [0.1; x0]));
plot(spec,[1 2],'--');

% formatting
xlim([0,0.1]);
ylim([0,0.06]);
xlabel('t');
ylabel('x');

text = ['Brake,BRKNC01,',num2str(res),',',num2str(tVer+tComp)];

% ------------------------------ END OF CODE ------------------------------
