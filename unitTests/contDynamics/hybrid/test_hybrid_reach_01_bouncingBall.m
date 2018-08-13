function res = test_hybrid_reach_01_bouncingBall
% test_hybrid_reach_01_bouncingBall - unit_test_function for hybrid
% dynamics
%
% Checks the solution of the hybrid system class for the bouncing ball
% example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been 
% saved. It is also checked whether the simulation matches the analytical
% solution.
%
% Syntax:  
%    res = test_hybrid_reach_01_bouncingBall
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff
% Written:      27-July-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


%set options---------------------------------------------------------------
options.x0 = [1; 0]; %initial state for simulation
options.R0 = zonotope([options.x0, diag([0.05, 0.05])]); %initial state for reachability analysis
options.startLoc = 1; %initial location
options.finalLoc = 0; %0: no final location
options.tStart = 0; %start time
options.tFinal = 1.7; %final time
options.timeStepLoc{1} = 0.05; %time step size for reachable set computation in location 1
options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.polytopeOrder = 10;
options.errorOrder=2;
options.reductionTechnique = 'girard';
options.guardIntersect = 'polytope';
options.isHyperplaneMap = 0;
options.enclosureEnables = 5; %choose enclosure method(s)
options.originContained = 0;
options.polytopeType = 'mpt';
%--------------------------------------------------------------------------


%specify hybrid automaton--------------------------------------------------
%specify linear system of bouncing ball
A = [0 1; 0 0];
B = [0; 0]; % no inputs
c = [0; -9.81]; % constant part of affine dynamics
linSys = linearSys('linearSys',A,B,c);

%define large and small distance
dist = 1e3;
eps = 1e-6;
alpha = -0.75; %rebound factor

%invariant
inv = interval([-2*eps; -dist], [dist; dist]);
%guard sets
guard = interval([-eps; -dist], [0; -eps]); 
%resets
reset.A = [0, 0; 0, alpha]; reset.b = zeros(2,1);
%transitions
trans{1} = transition(guard,reset,1,'a','b'); %--> next loc: 1; 'a', 'b' are dummies
%specify location
loc{1} = location('loc1',1,inv,trans,linSys); 
%specify hybrid automata
HA = hybridAutomaton(loc); % for "geometric intersection"
%--------------------------------------------------------------------------

%set input:
options.uLoc{1} = 0; % no inputs
options.uLocTrans{1} = 0; % no inputs
options.Uloc{1} = zonotope(0); % no inputs

%simulate hybrid automaton
HA = simulate(HA,options); 

%compute reachable set
[HA] = reach(HA,options);

%extract hitting times and states; remove artificial breaks due to guard
%detection
val = get(HA,'trajectory');
hitCounter = 0;
for iLoc = 1:length(val.x)-1
    if abs(val.x{iLoc}(end,1)) < 1e-10
        hitCounter = hitCounter + 1;
        tHit_sim(hitCounter) = val.t{iLoc+1}(1,:);
        xHit_sim(hitCounter,:) = val.x{iLoc+1}(1,:);
    end
end

%CHECK SIMULATION----------------------------------------------------------
%obtain exact hitting times and states
s0 = options.x0(1);
v0 = options.x0(2);
g = c(2);
tHit(1) = (-v0-sqrt(v0^2-2*g*s0))/g;
vHit(1) = alpha*(v0 + g*tHit(1));
xHit(1,:) = [0, vHit(1)];
for iLoc = 2:hitCounter
    tHit(iLoc) = tHit(iLoc-1) - 2/g*vHit(iLoc-1);
    vHit(iLoc) = -alpha*vHit(iLoc-1);
    xHit(iLoc,:) = [0, vHit(iLoc)];
end

%compare exact values with simulation results
res_tHit = all(abs(tHit - tHit_sim) < 1e-10);
res_xHit = all(all(abs(xHit - xHit_sim) < 1e-10));
res_sim = res_tHit & res_xHit;
%--------------------------------------------------------------------------

%CHECK REACHABLE SET-------------------------------------------------------
R = get(HA,'continuousReachableSet');
I = interval(R.OT{end}{end});

%saved result
I_saved = interval( ...
           [0.1481306674226159; -0.2059360624105376], ...
           [0.1974097982443197; 0.4812208662215985]);
        
%check if slightly bloated versions enclose each other
res_1 = (I <= enlarge(I_saved,1+1e-8));
res_2 = (I_saved <= enlarge(I,1+1e-8));

%final result
res_reach = res_1*res_2;
%--------------------------------------------------------------------------

%final result
res = res_sim & res_reach;

%------------- END OF CODE --------------
