function res = test_hybridAutomaton_reach_01_bouncingBall
% test_hybridAutomaton_reach_01_bouncingBall - unit test function for
%    hybrid systems with linear dynamics
%
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been 
% saved. It is also checked whether the simulation matches the analytical
% solution.
%
% Syntax:
%    res = test_hybridAutomaton_reach_01_bouncingBall
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff
% Written:       27-July-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Reachability Options ----------------------------------------------------

% problem description
params.R0 = zonotope([[1;0],diag([0.05,0.05])]);    % initial set
params.startLoc = 1;                                % initial location
params.finalLoc = 0;                                % 0: no final location
params.tFinal = 1.7;                                % final time

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'pca'};
options.intersectInvariant = true;

% Hybrid Automaton --------------------------------------------------------

% continuous dynamics 
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];
linSys = linearSys('linearSys',A,B,c);

% system parameters
alpha = -0.75;                                      % rebound factor

% invariant set 
inv = polytope([-1,0],0);

% guard sets
guard = conHyperplane([1,0],0,[0,1],0);

% reset function
reset.A = [0, 0; 0, alpha]; 
reset.c = zeros(2,1);

% transitions
trans = transition(guard,reset,1);

% location object
loc = location('loc1',inv,trans,linSys); 

% hybrid automaton
HA = hybridAutomaton(loc);


% Simulation and Reachability Analysis ------------------------------------

% compute reachable set
R = reach(HA,params,options);

% simulate hybrid automaton
params.x0 = center(params.R0);
params = rmfield(params,'R0');

[t,x] = simulate(HA,params); 
 
% extract hitting times and states; remove artificial breaks due to guard
% detection
tHit_sim = cellfun(@(x) x(1),t(2:6),'UniformOutput',true);
xHit_sim = cell2mat(cellfun(@(y) y(1,:),x(2:6),'UniformOutput',false));
tHit_sim = tHit_sim(xHit_sim(:,1) < 0.1);
xHit_sim = xHit_sim(xHit_sim(:,1) < 0.1,:);
hitCounter = length(tHit_sim);


% Check Simulation --------------------------------------------------------

% obtain exact hitting times and states from analytical solution
s0 = params.x0(1);
v0 = params.x0(2);
g = c(2);
tHit(1) = (-v0-sqrt(v0^2-2*g*s0))/g;
vHit(1) = alpha*(v0 + g*tHit(1));
xHit(1,:) = [0, vHit(1)];

for iLoc = 2:hitCounter
    tHit(iLoc,1) = tHit(iLoc-1) - 2/g*vHit(iLoc-1);
    vHit(iLoc) = -alpha*vHit(iLoc-1);
    xHit(iLoc,:) = [0, vHit(iLoc)];
end

% compare exact values with simulation results
res_tHit = all(withinTol(tHit,tHit_sim,1e-10));
res_xHit = all(all(withinTol(xHit,xHit_sim,1e-10)));
res_sim = res_tHit && res_xHit;


% Check Reachable Set -----------------------------------------------------

I = interval(R(end).timeInterval.set{end});

% saved result
I_saved = interval( ...
           [0.1481306674226159; -0.2059360624105376], ...
           [0.1974097982443197; 0.4812208662215985]);

% final result
res_reach = isequal(I,I_saved,1e-4);

% overall result
res = res_sim && res_reach;

% ------------------------------ END OF CODE ------------------------------
