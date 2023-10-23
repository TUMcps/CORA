function res = example_hybrid_reach_06_bouncingBallSineWave
% example_hybrid_reach_06_bouncingBallSineWave - example in Sec. 7.3.3 in
%    [1] describing a ball bouncing on a sine-wave shaped survace. The
%    resulting hybrid system has nonlienar guard set and nonlinear reset
%    functions
%
% Syntax:
%    res = example_hybrid_reach_06_bouncingBallSineWave
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] D. Ishii and et al. "An interval-based SAT modulo ODE solver for 
%        model checking nonlinear hybrid systems", 2011

% Authors:       Niklas Kochdumper
% Written:       10-December-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

R0 = interval([0.48;4.98;0;-5],[0.52;5.02;0;-5]);

params.R0 = zonotope(R0);      
params.startLoc = 1;                               
params.tFinal = 2.2;                 


% Reachability Options ----------------------------------------------------

% settings for continuous reachability 
options.alg = 'lin';
options.tensorOrder = 2;
options.timeStep = 0.01;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% settings for hybrid systems
options.guardIntersect = 'levelSet';


% Hybrid Automaton --------------------------------------------------------

% continuous dynamics 
f = @(x,u) [x(3); x(4); 0; -9.81+0.01*x(4)^2];
sys = nonlinearSys(f);

% invariant set 
syms x y vx vy
eq = -y + sin(x);
inv = levelSet(eq,[x;y;vx;vy],'<=');

% guard sets
guard = levelSet(-eq,[x;y;vx;vy],'==');

% reset function
reset.f = @(x,u) [x(1); ...
         x(2); ... 
         ((1-0.8*cos(x(1))^2)*x(3)+1.8*cos(x(1))*x(4))/(1+cos(x(1))^2); ...
         (1.8*cos(x(1))*x(3)+(-0.8+cos(x(1))^2)*x(4))/(1+cos(x(1))^2)];

% transitions
trans = transition(guard,reset,1);

% location object
loc = location(inv,trans,sys); 

% hybrid automaton
HA = hybridAutomaton(loc);


% Reachability Analysis ---------------------------------------------------

tic;
R = reach(HA,params,options);
tComp = toc;

disp(['Computation time for reachable set: ',num2str(tComp),' s']);


% Simulation --------------------------------------------------------------

simRes = simulateRandom(HA,params);


% Visualization -----------------------------------------------------------

figure; hold on;
xlim([-4,4]); ylim([-1,6]);
% axis([0,1.2,-6,4]);

% plot reachable set
plot(R);

% plot initial set
plot(params.R0,[1,2],'k','FaceColor','w');

% plot the ground
ground = levelSet(y-sin(x),[x;y;vx;vy],'<=');
plot(ground,[1,2],'FaceColor',colorblind('gray'));

% plot simulated trajectories
plot(simRes);

% completed successfully
res = true;

% ------------------------------ END OF CODE ------------------------------
