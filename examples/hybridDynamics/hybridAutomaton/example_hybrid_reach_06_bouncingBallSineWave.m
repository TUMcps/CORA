function res = example_hybrid_reach_06_bouncingBallSineWave
% example_hybrid_reach_01_bouncingBallSineWave - example in Sec. 7.3.3 in
%       [1] describing a ball bouncing on a sine-wave shaped survace. The
%       resulting hybrid system has nonlienar guard set and nonlinear reset
%       functions
%
% Syntax:  
%    example_hybrid_reach_06_bouncingBallSineWave
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean  
%
% References:
%    [1] D. Ishii and et al. "An interval-based SAT modulo ODE solver for 
%        model checking nonlinear hybrid systems", 2011

% Author:       Niklas Kochdumper
% Written:      10-December-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% Parameter ---------------------------------------------------------------

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
reset.f = @(x) [x(1);x(2); ... 
         ((1-0.8*cos(x(1))^2)*x(3)+1.8*cos(x(1))*x(4))/(1+cos(x(1))^2); ...
         (1.8*cos(x(1))*x(3)+(-0.8+cos(x(1))^2)*x(4))/(1+cos(x(1))^2)];

% transitions
trans{1} = transition(guard,reset,1);

% location object
loc{1} = location(inv,trans,sys); 

% hybrid automata
HA = hybridAutomaton(loc);


% Reachability Analysis ---------------------------------------------------

tic;
R = reach(HA,params,options);
tComp = toc;

disp(['Computation time for reachable set: ',num2str(tComp),' s']);


% Simulation --------------------------------------------------------------

% settings for random simulation
simOpt.points = 10;        % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.inpChanges = 10;    % changes of input over time horizon  

% random simulation
simRes = simulateRandom(HA,params,simOpt); 


% Visualization -----------------------------------------------------------

figure; hold on; box on; xlim([-4,4]); ylim([-1,6]);

% plot reachable set
plot(R,[1,2],'FaceColor',[0.7 0.7 0.7],'EdgeColor','none');

% plot initial set
plot(params.R0,[1,2],'FaceColor','w','EdgeColor','k','Filled',true);

% plot the ground
ground = levelSet(y-sin(x),[x;y;vx;vy],'<=');
plot(ground,[1,2],'FaceColor',[0 0 0.8],'Filled',true);

% plot simulated trajectories
plot(simRes,[1,2],'k');


res = 1;


%------------- END OF CODE --------------
