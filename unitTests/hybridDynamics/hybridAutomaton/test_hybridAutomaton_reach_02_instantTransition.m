function res = test_hybridAutomaton_reach_02_instantTransition()
% test_hybridAutomaton_reach_02_instantTransition - test for reachability
%    of hybrid dynamics, where the hybrid automaton contains instant
%    transitions (no elapsed time in between two subsequent transitions)
%
% Syntax:
%    test_hybridAutomaton_reach_02_instantTransition
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false  

% Authors:       Mark Wetzlinger
% Written:       17-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% Parameters --------------------------------------------------------------

% problem description
params.R0 = zonotope([0;-0.75],diag([0.05,0.05]));  % initial set
params.startLoc = 1;                                % initial location
params.tFinal = 10;                                 % final time


% Reachability Options ----------------------------------------------------

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 1;
options.zonotopeOrder = 5;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'box'};


% Hybrid Automaton --------------------------------------------------------

% simple construction: 2D automaton with three locations
% loc1: flow right
% loc2: flow upwards (irrelevant)
% loc3: flow left
% loc4: flow down (irrelevant)
% however, the transitions loc2->loc3 and loc4->loc1 are instant, thus
% we should only ever move from right to left and vice versa as the
% invariants in all locations are the same [-1,1] box

% continuous dynamics 
A = [0 0; 0 0];
B = [0; 0];
c1 = [1; 0];
c2 = [0; 1];
c3 = [-1; 0];
c4 = [0; -1];
linSys1 = linearSys('linearSys',A,B,c1);
linSys2 = linearSys('linearSys',A,B,c2);
linSys3 = linearSys('linearSys',A,B,c3);
linSys4 = linearSys('linearSys',A,B,c4);

% invariant set (same for all locations)
inv = polytope(interval([-1;-1],[1;1]));

% guard sets
guard1 = conHyperplane([1,0],1);
guard2 = fullspace(2);
guard3 = conHyperplane([-1,0],1);
guard4 = fullspace(2);

% reset functions
resetA = [1, 0; 0, 1];
resetc1 = [0; 0.15];
resetc2 = [-0.5; 0];
resetc3 = [0; 0.15];
resetc4 = [0.5; 0];
reset1 = struct('A',resetA,'c',resetc1);
reset2 = struct('A',resetA,'c',resetc2);
reset3 = struct('A',resetA,'c',resetc3);
reset4 = struct('A',resetA,'c',resetc4);

% transitions
trans1 = transition(guard1,reset1,2);
trans2 = transition(guard2,reset2,3);
trans3 = transition(guard3,reset3,4);
trans4 = transition(guard4,reset4,1);

% location objects
loc(1) = location('right',inv,trans1,linSys1);
loc(2) = location('up',inv,trans2,linSys2);
loc(3) = location('left',inv,trans3,linSys3);
loc(4) = location('down',inv,trans4,linSys4);

% instantiate hybrid automaton
HA1 = hybridAutomaton(loc);


% same hybrid automaton without instant transitions
% reset functions
resetc12 = [-0.5; 0.15];
resetc34 = [0.5; 0.15];
reset12 = struct('A',resetA,'c',resetc12);
reset34 = struct('A',resetA,'c',resetc34);

% transitions
trans1(1) = transition(guard1,reset12,2);
trans3(1) = transition(guard3,reset34,1);

% location objects
clear loc
loc(1) = location('right',inv,trans1,linSys1);
loc(2) = location('left',inv,trans3,linSys3);

HA2 = hybridAutomaton(loc);


% Reachability Analysis ---------------------------------------------------

R1 = reach(HA1,params,options);
R2 = reach(HA2,params,options);


% Numerical check ---------------------------------------------------------

% all sets need to be equal as well (only check first and last time point)
i_R1 = 1;
for i=1:length(R2)
    % skip index if from instant transition
    while isempty(R1(i_R1).timeInterval)
        i_R1 = i_R1 + 1;
    end

    % same number of sets per branch
    if length(R1(i_R1).timePoint.set) ~= length(R2(i).timePoint.set)
        res = false;
        break
    end

    % same start set and end set (faster than checking all, and should
    % catch any potential errors)
    if ~isequal(R1(i_R1).timePoint.set{1},R2(i).timePoint.set{1}) || ...
            ~isequal(R1(i_R1).timePoint.set{end},R2(i).timePoint.set{end})
        res = false;
        break
    end

    % increment index of R2
    i_R1 = i_R1 + 1;
end

% ------------------------------ END OF CODE ------------------------------
