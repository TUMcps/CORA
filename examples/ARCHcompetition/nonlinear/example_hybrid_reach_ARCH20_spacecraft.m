function res = example_hybrid_reach_ARCH20_spacecraft()
% example_hybrid_reach_ARCH20_spacecraft- spacecraft rendezvous benchmark
%                                         from the 2020 ARCH competition
%
% Syntax:  
%    res = example_hybrid_reach_ARCH20_spacecraft
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% 
% Author:       Niklas Kochdumper
% Written:      31-May-2018
% Last update:  13-March-2019
% Last revision:---


%------------- BEGIN CODE --------------

% x_1 = v_x
% x_2 = v_y
% x_3 = s_x 
% x_4 = s_y
% x_5 = t



% Parameter ---------------------------------------------------------------

% problem description
R0 = zonotope([[-900;-400;0;0;0],diag([25;25;0;0;0])]);

params.R0 = R0;                                     % initial set
params.startLoc = 1;                                % initial location
params.finalLoc = 0;                                % 0: no final location
params.tFinal = 200;                                % final time



% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.timeStep{1} = 2e-1;
options.timeStep{2} = 5e-2;
options.timeStep{3} = 2e-1;

options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 30;
options.zonotopeOrder = 10;

% settings for hybrid systems
options.guardIntersect = 'nondetGuard';
options.enclose = {'pca'}; 


% System Dynamics ---------------------------------------------------------

% specify hybrid automata (automatically converted from SpaceEx)
HA = rendeszvous_nonlinear_passive_nondet();


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

% passive mode -> check if the space station was hit
spaceStation = interval([-0.1;-0.1],[0.1;0.1]);
func = @(x) checkSpaceStationHit(x,spaceStation);

spec3 = specification(func,'custom');

% specifications for each mode
spec = cell(3,1);
spec{2} = add(spec1,spec2);
spec{3} = spec3;


% Reachability Analysis ---------------------------------------------------

% reachable set computations
tic
[R,res] = reach(HA,params,options,spec);
tComp = toc;

% display results
disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tComp)]);


% Visualiztion ------------------------------------------------------------

figure; hold on; box on

% plot line-of-sight cone
h = fill([-100,0,-100],[-60,0,60],'g','FaceAlpha',0.6,'EdgeColor','none');
set(h,'FaceColor','y');

% plot reachable set
R1 = find(R,'location',1);
plot(R1,[1,2],'b','Filled',true,'EdgeColor','none');

R1 = find(R,'location',2);
plot(R1,[1,2],'m','Filled',true,'EdgeColor','none');

R1 = find(R,'location',3);
plot(R1,[1,2],'g','Filled',true,'EdgeColor','none');

% plot initial set
plot(params.R0,[1,2],'w','Filled',true,'EdgeColor','k');

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