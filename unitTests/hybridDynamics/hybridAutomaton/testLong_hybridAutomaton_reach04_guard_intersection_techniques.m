function res = testLong_hybridAutomaton_reach04_guard_intersection_techniques
% testLong_hybridAutomaton_reach04_guard_intersection_techniques - unit test function for
%    hybrid systems with linear and nonlinear dynamics
%
% It is checked whether different guard intersection algorithms terminate.
% A future test will also check the quality of the results.
% We require both linear and nonlinear systems to be checked in this
% regard, since the guard intersections call the respective reach-functions
% internally.
%
% Syntax:
%    res = testLong_hybridAutomaton_reach04_guard_intersection_techniques
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Maximilian Perschl
% Written:       26-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;


% Reachability Options ----------------------------------------------------

% problem description
params.R0 = zonotope([[1;0],diag([0.05,0.05])]);    % initial set
params.startLoc = 1;                                % initial location
params.finalLoc = 0;                                % 0: no final location
params.tFinal = 1;                                  % final time

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% settings for hybrid systems
% options.guardIntersect = 'pancake';
% options.guardOrder = 3;
% options.intersectInvariant = true;
options.alg = 'lin';
options.tensorOrder = 2;

% Hybrid Automaton --------------------------------------------------------

% continuous dynamics 
A = [0 1; 0 0];
B = [0; 0];
c = [0; -9.81];

linsys = linearSys('linearSys',A,B,c);
nonlinSys = nonlinearSys(linsys);

% system parameters
alpha = -0.75;                                      % rebound factor

% invariant set 
inv = polytope([-1,0],0);

% guard sets
guard = polytope([0,1],0,[1,0],0);

% reset function
reset = linearReset([0, 0; 0, alpha], zeros(2,1), zeros(2,1)); 

% transitions
trans = transition(guard,reset,1);

% location object
loc_lin = location('loc1',inv,trans,linsys); 
loc_nonlin = location('loc1',inv,trans,nonlinSys); 

% hybrid automaton
HA_lin = hybridAutomaton(loc_lin);
HA_nonlin = hybridAutomaton(loc_nonlin);

% Reachability Analysis ---------------------------------------------------

guardIntersect = {'polytope','zonoGirard','conZonotope', ...
                  'hyperplaneMap','pancake','nondetGuard'};

for i = 1:length(guardIntersect)

    % set options for guard intersection 
    options_ = options;
    options_.guardIntersect = guardIntersect{i};
    
    switch guardIntersect{i}
        % execute different intersection algorithms/parameters
        case 'polytope'
            options_.enclose = {'pca'};
            
        case 'zonoGirard'
            options_.enclose = {'pca'};
            
        case 'conZonotope'
            options_.enclose = {'pca'};
            options_.guardOrder = 5;
            
        case 'hyperplaneMap'
            options_.guardOrder = 3;
            
        case 'nondetGuard'
            options_.enclose = {'pca'};
    end

    % compute reachable set
    R_lin = reach(HA_lin,params,options_);
    assertLoop(isa(R_lin,"reachSet") && query(R_lin,"tFinal") > 0,i);
    R_nonlin = reach(HA_nonlin,params,options_);
    assertLoop(isa(R_nonlin,"reachSet") && query(R_nonlin,"tFinal") > 0,i);
end

% ------------------------------ END OF CODE ------------------------------
