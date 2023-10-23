function res = test_parallelHybridAutomaton_reach_01_emptyInvGuard()
% test_parallelHybridAutomaton_reach_01_emptyInvGuard - test for
%    reachability of hybrid dynamics, where the parallel hybrid automaton
%    contains empty invariants and guard sets
%
% Syntax:
%    test_parallelHybridAutomaton_reach_01_emptyInvGuard
%
% Inputs:
%    no
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       25-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

% parallel hybrid automaton
params.R0 = zonotope([0;-0.75;0.75;0],diag(0.05*ones(4,1)));
params.U = zonotope(0);
params.startLoc = [1;1];
params.tFinal = 10;

% hybrid automaton
paramsHA.R0 = zonotope([0;-0.75],diag(0.05*ones(2,1)));
paramsHA.U = zonotope(0);
paramsHA.startLoc = 1;
paramsHA.tFinal = 10;

% Reachability Options ----------------------------------------------------

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 1;
options.zonotopeOrder = 5;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'box'};


% Parallel Hybrid Automaton -----------------------------------------------

% idea: two hybrid automata in 2D
% invariant is a box, flows are linear
% HA1: flow right and left
% HA2: flow up and down (but slower than HA1)
% -> whenever HA1 hits right or left, HA2 switches location as well
%    via instant transition/synchronization label

% continuous dynamics: right/left faster than up/down
A = [0 0; 0 0];
B = [0; 0];
c_right = [1; 0];
c_up = [0; 0.5];
c_left = [-1; 0];
c_down = [0; -0.5];
linSys_right = linearSys('right',A,B,c_right);
linSys_up = linearSys('up',A,B,c_up);
linSys_left = linearSys('left',A,B,c_left);
linSys_down = linearSys('down',A,B,c_down);

% invariant set 
inv = polytope(interval([-1;-1],[1;1]));
inv_fullspace = fullspace(2);

% guard sets (boundary of box)
guard_right = conHyperplane([1,0],1);
guard_left = conHyperplane([-1,0],1);
% guard set for instant transition
guard_instant = fullspace(2);

% reset functions: shift up/left
resetA = [1, 0; 0, 1];
resetc_up = [0; 0.15];
resetc_left = [-0.15; 0];
reset_up = struct('A',resetA,'c',resetc_up);
reset_left = struct('A',resetA,'c',resetc_left);

% HA1: transitions
trans1 = transition(guard_right,reset_up,2,'hit_right');
trans2 = transition(guard_left,reset_up,1,'hit_left');

% HA1: locations
loc = [location('right',inv,trans1,linSys_right);...
        location('left',inv,trans2,linSys_left)];

% instantiate first hybrid automaton
HA1 = hybridAutomaton(loc);

% HA2: transitions
trans1 = transition(guard_instant,reset_left,2,'hit_right');
trans2 = transition(guard_instant,reset_left,1,'hit_left');

% HA2: locations (invariants are full space)
loc = [location('up',inv_fullspace,trans1,linSys_up);...
        location('down',inv_fullspace,trans2,linSys_down)];

% instantiate second hybrid automaton
HA2 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2];
inputBinds{1} = [0 1];
inputBinds{2} = [0 1];

% instantiate parallel hybrid automaton
pHA = parallelHybridAutomaton(components,inputBinds);


% Reachability Analysis ---------------------------------------------------

% parallel hybrid automaton
R = reach(pHA,params,options);

% first hybrid automaton
R_HA = reach(HA1,paramsHA,options);


% Visualization -----------------------------------------------------------
% just for debugging purposes

% figure;
% projDims = {[1,2],[3,4]};
% for p=1:2
%     subplot(1,2,p); hold on; box on;
% 
%     % invariant
%     if p == 1
%         plot(inv,[1,2],'k');
%     end
% 
%     % reachable set
%     plot(R,projDims{p});
%     % reachable set (only first projection)
%     if p == 1
%         plot(R_HA,[1,2],'EdgeColor',colorblind('r'));
%     end
% 
%     % initial set
%     plot(project(params.R0,projDims{p}),[1,2],...
%         'FaceColor','w','EdgeColor','k');
% 
%     % axes and labels
%     axis([-1.25,1.25,-1.25,1.25]);
%     xlabel("x_" + projDims{p}(1));
%     ylabel("x_" + projDims{p}(2));
% end


% Numerical check ---------------------------------------------------------

% since instant transitions do not effect additional branches in the
% reachSet object, both reachSet objects have to have the same length
res = length(R) == length(R_HA);

% ensure that start sets after transitions are equal
for i=1:length(R)
    if ~isequal(project(R(i).timePoint.set{1},[1,2]),...
            R_HA(i).timePoint.set{1},1e-14)
        res = false;
        break
    end
end

% ------------------------------ END OF CODE ------------------------------
