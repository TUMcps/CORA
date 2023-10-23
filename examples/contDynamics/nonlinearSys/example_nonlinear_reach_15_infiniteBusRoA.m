function completed = example_nonlinear_reach_15_infiniteBusRoA
% example_nonlinear_reach_15_infiniteBusRoA - example for verifying
%    the region of attraction of a single-machine-infinite-bus system from
%    [1]
%
% Syntax:
%    completed = example_nonlinear_reach_15_infiniteBusRoA
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% References:
%    [1] M. Althoff. "Formal Verification of Power Systems", submitted to
%        ARCH 2022.

% Authors:       Matthias Althoff
% Written:       02-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 2.9;
params.R0 = zonotope([[asin(1/5); 0],0.7*eye(2)]);
params.U = zonotope([0, 0]);


% Reachability Settings ---------------------------------------------------

% up to 0.8:
options.timeStep = 0.002;
options.taylorTerms = 4;
options.zonotopeOrder = 200;
options.alg = 'lin';
options.tensorOrder = 3;
options.errorOrder = 2;
options.intermediateOrder = 2;


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@infiniteBus);


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Check Enclosure in Invariant --------------------------------------------

% rotation angle and matrix
%rotAngle = 0.003;
rotAngle = 0.005;
T = [cos(rotAngle) -sin(rotAngle); sin(rotAngle) cos(rotAngle)];

% unrotated ellipsoid
Q = diag([0.02,5]);
Eorig = ellipsoid(Q,[asin(1/5); 0]);

%rotated ellipsoid
E = T*Eorig;
enclosure = contains(E, R.timePoint.set{end});
disp("Enclosure satisfied? " + all(enclosure));

%% Simulation 

simOpt.points = 60;
%simOpt.points = 6;
simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

figure; hold on; box on;
useCORAcolors("CORA:contDynamics")

% plot reachable sets
plot(R);

% plot initial set
plot(R(1).R0,[1 2]);

% plot simulation results     
plot(simRes);

% plot invariant
plot(E,[1 2],'LineWidth',3,'Color',colorblind('g'));

% label plot
xlabel('a');
ylabel('b');


% plot final reachable set
figure; hold on; box on;

plot(specification(E, 'safeSet'),[1 2]);
plot(R.timePoint.set{end},[1,2],'FaceColor',CORAcolor("CORA:reachSet"));


% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
