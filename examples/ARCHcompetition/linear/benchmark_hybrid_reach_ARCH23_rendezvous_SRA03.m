function text = benchmark_hybrid_reach_ARCH23_rendezvous_SRA03()
% benchmark_hybrid_reach_ARCH23_rendezvous_SRA03 - spacecraft rendezvous
%    example from the 2023 ARCH competition
%
% Syntax:
%    res = benchmark_hybrid_reach_ARCH23_rendezvous_SRA03()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission"  

% Authors:       Niklas Kochdumper
% Written:       20-March-2018
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% x_1 = v_x
% x_2 = v_y
% x_3 = s_x 
% x_4 = s_y
% x_5 = t


% Parameter ---------------------------------------------------------------

R0 = zonotope([[-900; -400; 0; 0; 0],diag([25; 25; 0; 0; 0])]);
params.R0 = R0;

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
options.enclose = {'pca'};


% System Dynamics ---------------------------------------------------------

HA = rendezvous_SRA03();


% Specification -----------------------------------------------------------

% randezvous attempt -> check if velocity inside feasible region
C = [0 0 1 0 0;0 0 -1 0 0;0 0 0 1 0;0 0 0 -1 0; ....
     0 0 1 1 0;0 0 1 -1 0;0 0 -1 1 0;0 0 -1 -1 0];
d = [3;3;3;3;4.2;4.2;4.2;4.2];
poly = polytope(C,d);
active_locs = 2;

spec1 = specification(poly,'safeSet',active_locs);

% randezvous attempt -> check if spacecraft inside line-of-sight
C = [tan(pi/6) -1 0 0 0;tan(pi/6) 1 0 0 0];
d = [0;0];
poly = polytope(C,d);
active_locs = 2;

spec2 = specification(poly,'safeSet',active_locs);

% passive mode -> check if the space station was hit
spaceStation = interval([-0.1;-0.1],[0.1;0.1]);
func = @(x) aux_checkSpaceStationHit(x,spaceStation);
active_locs = 3;

spec3 = specification(func,'custom',active_locs);

% specifications
spec = add(spec1,spec2);
spec = add(spec,spec3);


% Simulation --------------------------------------------------------------

simRes = simulateRandom(HA,params); 


% Reachability Analysis ---------------------------------------------------

% reachable set computations
timer = tic;
[R,res] = reach(HA,params,options,spec);
tComp = toc(timer);

% display results
disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tComp)]);


% Visualization -----------------------------------------------------------

figure; hold on; box on

% plot line-of-sight cone
h = fill([-100,0,-100],[-60,0,60],'g','EdgeColor','k');
set(h,'FaceColor',colorblind('gray'));

useCORAcolors('CORA:contDynamics')

% plot reachable set
plot(R,[1,2]);

% plot initial set
plot(R(1).R0,[1,2]);

% plot simulation
plot(simRes,[1,2]);

% plot space station
plot(specification(spaceStation),[1,2]);

xlabel('s_x');
ylabel('s_y');

text = ['Rendezvous,SRA03,',num2str(res),',',num2str(tComp)];

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkSpaceStationHit(S,spaceStation)
% check if the reachable set intersects the space station

   temp = project(S,[1,2]);
   res = true;
   
   % first check: is interval enclosure intersecting?
   if isIntersecting(spaceStation,interval(temp))
       
       % second check: is reduced zonotope intersecting?
       temp = reduce(temp,'girard',3);
       if isIntersecting(spaceStation,temp)
           res = false;
       end
   end
end

% ------------------------------ END OF CODE ------------------------------
