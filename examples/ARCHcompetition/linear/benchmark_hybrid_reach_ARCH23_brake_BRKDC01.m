function text = benchmark_hybrid_reach_ARCH23_brake_BRKDC01()
% benchmark_hybrid_reach_ARCH23_brake_BRKDC01 -  brake benchmark from the
%    2023 ARCH competition
%
% Syntax:
%    benchmark_hybrid_reach_ARCH23_brake_BRKDC01
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

% system matrices
A = [-(R+K^2/d_rot)/L, 0, K_P/L, K_I/L; ....
     K/(i*d_rot), 0, 0, 0; ...
     0, 0, 0, 0; ...
     0, 0, 0, 0];
 
B = [0;0;0;0];

% reset
reset.A = [1 0 0 0; 0 1 0 0; 0 -1 0 0;0 -T_sample 0 1];
reset.c = [0; 0; x0; T_sample*x0];

% linear system object
sys = linearSys(A,B);


% Parameter ---------------------------------------------------------------

params.tFinal = T_sample;
params.R0 = zonotope([0;0;0;0]);


% Settings ----------------------------------------------------------------

options.timeStep = T_sample/5;
options.taylorTerms = 5;
options.zonotopeOrder = 20;


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = {};

for i = 1:1000
    
   % compute reachable set
   Rtemp = reach(sys,params,options);
   R = [R; Rtemp.timeInterval.set];
   
   % compute reset
   params.R0 = reset.A*Rtemp.timePoint.set{end} + reset.c;
    
end

tComp = toc(clock);

% instantiate reachSet object
t = 0;
timeInt.set = cell(length(R),1);
timeInt.time = cell(length(R),1);
for i = 1:length(R)
    
    % get intervals
    timeInt.set{i} = interval(project(R{i},2));
    timeInt.time{i} = interval(t,t + options.timeStep);
    
    % update time
    t = t + options.timeStep;
end

R_set = reachSet([],timeInt);


% Verification ------------------------------------------------------------

clock = tic;
res = true;

for i = 1:length(R)
   temp = supremum(interval(project(R{i},2)));
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

useCORAcolors("CORA:contDynamics");

% plot reachable set
plotOverTime(R_set,1,'Unify',true);

% plot specification
spec = specification(interval([0; x0], [0.1; x0]));
plot(spec,[1 2],'--');

% formatting
xlim([0,0.1]);
ylim([0,0.06]);
xlabel('$t$','interpreter','latex');
ylabel('$x$','interpreter','latex');


text = ['Brake,BRKDC01,',num2str(res),',',num2str(tVer+tComp)];

% ------------------------------ END OF CODE ------------------------------
