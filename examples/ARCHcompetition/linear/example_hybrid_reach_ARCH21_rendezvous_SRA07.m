function res = example_hybrid_reach_ARCH21_rendezvous_SRA07()
% example_hybrid_reach_ARCH21_rendezvous_SRA07 - spacecraft rendezvous
%    benchmark from the ARCH 2021 competition
%
% Syntax:  
%    res = example_hybrid_reach_ARCH21_rendezvous_SRA07()
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
% Written:      12-March-2018
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
options.timeStep{3} = 2e-1;

% settings for continuous reachability 
options.zonotopeOrder=10; 
options.taylorTerms=3;

% settings for computing guard intersections
options.guardIntersect = 'nondetGuard';
options.enclose = {'pca','box'};


% System Dynamics ---------------------------------------------------------

HA = rendezvous_SRA07();


% Specification -----------------------------------------------------------

% randezvous attempt -> check if velocity inside feasible region
C = [0 0 1 0 0;0 0 -1 0 0;0 0 0 1 0;0 0 0 -1 0; ....
     0 0 1 1 0;0 0 1 -1 0;0 0 -1 1 0;0 0 -1 -1 0];
d = [3;3;3;3;4.2;4.2;4.2;4.2];

func = @(x) ~any(supremum(interval(C*x) - d) > 0);

spec1 = specification(func,'custom');

% randezvous attempt -> check if spacecraft inside line-of-sight
C = [tan(pi/6) -1 0 0 0;tan(pi/6) 1 0 0 0];
d = [0;0];

func = @(x) ~any(supremum(interval(C*x) - d) > 0);

spec2 = specification(func,'custom');

% passive mode -> check if the space station was hit
spaceStation = interval([-0.1;-0.1],[0.1;0.1]);
func = @(x) checkSpaceStationHit(x,spaceStation);

spec3 = specification(func,'custom');

% specifications for each mode
spec = cell(3,1);
spec{2} = add(spec1,spec2);
spec{3} = spec3;


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

% plot space station
plot(spaceStation,[1,2],'r','Filled',true,'EdgeColor','none');

xlabel('s_x');
ylabel('s_y');

end


% Auxiliary Functions -----------------------------------------------------

function res = checkSpaceStationHit(S,spaceStation)
% check if the reachable set intersects the space station

   temp = project(S,[1,2]);
   res = 1;
   
   % first check: is interval enclosure intersecting?
   if isIntersecting(spaceStation,interval(temp))
       
       % second check: is reduced zonotope intersecting?
       temp = reduce(temp,'girard',3);
       if isIntersecting(spaceStation,temp)
           res = 0;
       end
   end
end

%------------- END OF CODE --------------