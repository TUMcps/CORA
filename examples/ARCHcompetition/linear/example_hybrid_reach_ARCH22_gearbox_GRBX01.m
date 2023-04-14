function example_hybrid_reach_ARCH22_gearbox_GRBX01()
% example_hybrid_reach_ARCH22_gearbox_GRBX01 - gearbox benchmark from the
%                                              2022 ARCH competition
%
% Syntax:  
%    example_hybrid_reach_ARCH22_gearbox_GRBX01
%
% Inputs:
%    -
%
% Outputs:
%    -

% Author:       Niklas Kochdumper
% Written:      29-May-2018
% Last update:  13-March-2019
% Last revision:---

%------------- BEGIN CODE --------------


% Parameters --------------------------------------------------------------

% initial set
int = interval([0;0;-0.0168;0.0029;0],[0;0;-0.0166;0.0031;0]);
params.R0 = zonotope(int);

% initial location, final location, and final time
params.startLoc = 1;
params.finalLoc = 2;
params.tFinal = 0.2;


% Reachability Options ----------------------------------------------------

% settings for continuous reachability
options.timeStep = 0.2/182; %1.1e-3 does not divide tFinal in integer steps
options.zonotopeOrder = 20;
options.taylorTerms = 3;

% settings for computing guard intersections
options.guardIntersect = 'zonoGirard';
options.enclose = {'box','pca'};


% Hybrid Automaton --------------------------------------------------------

HA = gearbox();


% Reachability Analysis ---------------------------------------------------

clock = tic;
R = reach(HA,params,options);
tComp = toc(clock);


% Verification ------------------------------------------------------------

list = query(R,'reachSet');
clock = tic;

% constraint x5 < 20
res = 1;

for i = 1:length(list)
    temp = interval(project(list{i},5));
    if supremum(temp) >= 20
        res = 0;
        break;
    end
end

tVer = toc(clock);

disp(['specifications verified: ',num2str(res)]);
disp(['computation time: ',num2str(tVer + tComp)]);

% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot reachable set
useCORAcolors("CORA:contDynamics")
plot(R,[3,4]);
plot(R(1).R0,[3,4])

% plot tooth of gear
locs = HA.location;
trans = locs{1}.transition;
guard1 = trans{1}.guard;
guard2 = trans{2}.guard;

specs = specification({ ...
    halfspace(guard1.h.c(3:4),guard1.h.d); ...
    halfspace(guard2.h.c(3:4),guard2.h.d) ...
}, 'unsafeSet');

axis([-2e-2,-3e-3,-1e-2,1e-2]);

plot(specs,[1,2]);

%------------- END OF CODE --------------