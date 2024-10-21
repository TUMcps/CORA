function res = test_hybridAutomaton_reach_unsafeSet
% test_hybridAutomaton_reach_unsafeSet - test function for reachability
%    analysis of hybrid systems; here, we ensure that a reachable set that
%    has already left the invariant cannot intersect an unsafe set
%
% Syntax:
%    res = test_hybridAutomaton_reach_unsafeSet
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       21-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% simple automaton: the reachable set moves from left to right and vice
% versa between two vertical guards, behind which there are unsafe sets
inv = polytope(interval([-1;-100],[1;100]));
guard_right = polytope([],[],[1 0],1);
guard_left = polytope([],[],[1 0],-1);
reset = linearReset(eye(2),zeros(2,1),[0;1]);
trans1 = transition(guard_right,reset,2);
trans2 = transition(guard_left,reset,1);

% dynamics
linsys1 = linearSys(zeros(2,2),zeros(2,1),[1;0]);
linsys2 = linearSys(zeros(2,2),zeros(2,1),[-1;0]);

% init locations and hybrid automaton
loc1 = location('to-right',inv,trans1,linsys1);
loc2 = location('to-left',inv,trans2,linsys2);
HA = hybridAutomaton([loc1;loc2]);

% model parameters: make initial set wide enough so that it reaches the
% unsafe set behind the guard set before fully exiting the invariant set;
% also, make the time horizon long enough so that there are a couple of
% bounces left and right
params.R0 = zonotope([0;0],[0.75 0; 0 0.25]);
params.startLoc = 1;
params.finalLoc = 0;
params.tFinal = 5;

% settings for continuous reachability 
options.timeStep = 0.05;
options.taylorTerms = 10;
options.zonotopeOrder = 20;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'pca'};

% specifications: unsafe sets
S_right = interval([1.5;-0.1],[2;0.1]);
S_left = interval([-2;-0.1],[-1.5;0.1]);
spec = [specification(S_right,'unsafeSet',1); ...
        specification(S_left,'unsafeSet',2)];

% compute reachable set
R = reach(HA,params,options,spec);

% final set must reach time horizon (no abortion due to unsafe sets)
assert(contains(R(end).timePoint.time{end},params.tFinal));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
