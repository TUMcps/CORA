function res = example_hybrid_reach_ARCH21_rendezvous_SRNA01()
% example_hybrid_reach_ARCH21_rendezvous_SRNA01 - spacecraft rendezvous
%    benchmark from the ARCH 2021 competition
%
% Syntax:  
%    res = example_hybrid_reach_ARCH21_rendezvous_SRNA01()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission"  

% Author:       Niklas Kochdumper
% Written:      12-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% x_1 = v_x
% x_2 = v_y
% x_3 = s_x 
% x_4 = s_y
% x_5 = t


% Parameter ---------------------------------------------------------------

params.R0 = zonotope([[-900; -400; 0; 0; 0],diag([25; 25; 0; 0; 0])]);
params.startLoc = 1;
params.tFinal = 300;


% Reachability Settings ---------------------------------------------------

% time step
options.timeStep{1} = 2e-1;
options.timeStep{2} = 2e-2;

% settings for continuous reachability 
options.zonotopeOrder=10; 
options.taylorTerms=3;

% settings for computing guard intersections
options.guardIntersect = 'zonoGirard';
options.enclose = {'box'};


% System Dynamics ---------------------------------------------------------

HA = rendezvous_SRNA01();


% Specification -----------------------------------------------------------

% randezvous attempt -> check if velocity inside feasible region
C = [0 0 1 0 0;0 0 -1 0 0;0 0 0 1 0;0 0 0 -1 0; ....
     0 0 1 1 0;0 0 1 -1 0;0 0 -1 1 0;0 0 -1 -1 0];
d = [3;3;3;3;4.2;4.2;4.2;4.2];
poly = mptPolytope(C,d);

spec1 = specification(poly,'safeSet');

% randezvous attempt -> check if spacecraft inside line-of-sight
C = [tan(pi/6) -1 0 0 0;tan(pi/6) 1 0 0 0];
d = [0;0];
poly = mptPolytope(C,d);

spec2 = specification(poly,'safeSet');

% specifications for each mode
spec = cell(3,1);
spec{2} = add(spec1,spec2);


% Simulation --------------------------------------------------------------

% settings for simulation
simOpt.points = 10;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 7;

% random simulation
simRes = simulateRandom(HA,params,simOpt); 


% Reachability Analysis ---------------------------------------------------

% reachable set computations
timer = tic;
[R,res] = reach(HA,params,options,spec);
tComp = toc(timer);

% display results
disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tComp)]);


% Visualiztion ------------------------------------------------------------

figure; hold on; box on

% plot line-of-sight cone
h = fill([-100,0,-100],[-60,0,60],'g','FaceAlpha',0.6,'EdgeColor','none');
set(h,'FaceColor',[0 .6 0]);

% plot reachable set
plot(R,[1,2],'FaceColor',[.6 .6 .6],'Filled',true,'EdgeColor','none');

% plot initial set
plot(params.R0,[1,2],'w','Filled',true,'EdgeColor','k');

% plot simulation
plot(simRes,[1,2]);

xlabel('s_x');
ylabel('s_y');

%------------- END OF CODE --------------