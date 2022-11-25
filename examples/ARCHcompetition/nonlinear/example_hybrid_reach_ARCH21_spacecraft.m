function res = example_hybrid_reach_ARCH21_spacecraft()
% example_hybrid_reach_ARCH21_spacecraft- spacecraft rendezvous benchmark
%                                         from the 2021 ARCH competition
%
% Syntax:  
%    res = example_hybrid_reach_ARCH21_spacecraft
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
R0 = zonotope([[-900;-400;2.5;2.5;0],diag([25;25;2.5;2.5;0])]);

params.R0 = R0;                                     % initial set
params.startLoc = 1;                                % initial location
params.finalLoc = 0;                                % 0: no final location
params.tFinal = 200;                                % final time



% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.timeStep{1} = 0.2;
options.timeStep{2} = 0.05;
options.timeStep{3} = 0.2;

options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 10;
options.zonotopeOrder = 10;

% settings for hybrid systems
options.guardIntersect = 'nondetGuard';
options.enclose = {'pca'}; 


% System Dynamics ---------------------------------------------------------

% specify hybrid automata (automatically converted from SpaceEx)
HA = rendeszvous_nonlinear_passive_nondet();


% Specification -----------------------------------------------------------

% rendezvous attempt -> check if velocity inside feasible region
C = [0 0 1 0 0;0 0 -1 0 0;0 0 0 1 0;0 0 0 -1 0; ....
     0 0 1 1 0;0 0 1 -1 0;0 0 -1 1 0;0 0 -1 -1 0];
d = [3;3;3;3;4.2;4.2;4.2;4.2];
poly = mptPolytope(C,d);

spec1 = specification(poly,'safeSet');

% rendezvous attempt -> check if spacecraft inside line-of-sight
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


% Visualization -----------------------------------------------------------

figure; hold on; box on;
sxsy = [1,2];

% plot line-of-sight cone
h = fill([-100,0,-100],[-60,0,60],'g','FaceAlpha',0.6,'EdgeColor','none');
set(h,'FaceColor','y');

% plot reachable set
R1 = find(R,'location',1);
plot(R1,sxsy,'b','Filled',true,'EdgeColor','none','Unify',true);

R2 = find(R,'location',2);
plot(R2,sxsy,'m','Filled',true,'EdgeColor','none','Unify',true);

R3 = find(R,'location',3);
plot(R3,sxsy,'g','Filled',true,'EdgeColor','none','Unify',true);

% plot initial set
plot(params.R0,sxsy,'w','Filled',true,'EdgeColor','k');

% plot space station
plot(spaceStation,[1,2],'c','Filled',true,'EdgeColor','none');

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