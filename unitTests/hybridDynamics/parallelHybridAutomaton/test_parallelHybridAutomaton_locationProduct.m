function res = test_parallelHybridAutomaton_locationProduct
% test_parallelHybridAutomaton_locationProduct - test function for
%    computation of location product
%
% Syntax:
%    res = test_parallelHybridAutomaton_locationProduct
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
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% first component, first location
dynamics{1,1} = linearSys([0 -1; 1 0],[1 0; 0 1],[],[0.05 0.05]);
inv{1,1} = polytope([1 1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[3;3]);
trans = transition(guard,reset,2);
loc = location('loc1',inv{1,1},trans,dynamics{1,1});

% first component, second location
dynamics{1,2} = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
inv{1,2} = polytope([-1 -1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[-3;3]);
trans = transition(guard,reset,1);
loc(2) = location('loc2',inv{1,2},trans,dynamics{1,2});

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
dynamics{2,1} = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
inv{2,1} = polytope([1 1 1],1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[1;1;1]);
trans = transition(guard,reset,2);
loc(1) = location('loc1',inv{2,1},trans,dynamics{2,1});

% second component, second location
dynamics{2,2} = linearSys([0 -1 1; 1 0 0; 0 1 0],[0;0;1],[],[0.05 0 0; 0 0.05 -0.05]);
inv{2,2} = polytope([-1 -1 -1],-1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[-1;-1;-1]);
trans = transition(guard,reset,1);
loc(2) = location('loc2',inv{2,2},trans,dynamics{2,2});

% init second hybrid automaton
HA2 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2];
inputBinds{1} = [2 1; 2 2];
inputBinds{2} = [1 1];

% instantiate parallel hybrid automaton
pHA = parallelHybridAutomaton(components,inputBinds);


% compute location product
locID = [1;2];
allLabels = struct('name',[],'component',[],...
            'location',[],'transition',[],'guardempty',[]);

loc = locationProduct(pHA,locID,allLabels);

% resulting invariant dimension is 5
res = dim(loc.invariant) == 5;
% compare with true invariant
inv_true = polytope([1,1,0,0,0;0,0,-1,-1,-1],[0;-1]);
res(end+1,1) = loc.invariant == inv_true;
% resulting flow equation is linear
res(end+1,1) = isa(loc.contDynamics,'linearSys');
% compare with true linear system
dynamics_true = linearSys('linearSys x linearSys',...
    [[0 -1; 1 0], [0.05 0 0; 0 0.05 -0.05];...
    [0,0;0,0;0.05,0.05], [0 -1 1; 1 0 0; 0 1 0]],zeros(5,1),zeros(5,1));
res(end+1,1) = isequal(loc.contDynamics,dynamics_true);
% two transitions
res(end+1,1) = length(loc.transition) == 2;
% targets
res(end+1,1) = all(loc.transition(1).target == [2;2]) ...
    && all(loc.transition(2).target == [1;1]);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
