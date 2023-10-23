function res = test_parallelHybridAutomaton_isequal
% test_parallelHybridAutomaton_isequal - test function for isequal
%
% Syntax:
%    res = test_parallelHybridAutomaton_isequal
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

% empty parallel hybrid automata
res = isequal(parallelHybridAutomaton(),parallelHybridAutomaton());

% first component, first location
dynamics{1,1} = linearSys([0 -1; 1 0],[1 0; 0 1],[],[0.05 0.05]);
inv{1,1} = polytope([1 1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[3;3]);
trans = transition(guard,reset,2);
loc = location('loc1',inv{1,1},trans,dynamics{1,1});

% first component, second location
dynamics{1,2} = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
dynamics_ = linearSys([0 -1; -1 0],[0 1; 2 0],[],[0.05 -0.05]);
inv{1,2} = polytope([-1 -1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[-3;3]);
trans = transition(guard,reset,1);
loc(2) = location('loc2',inv{1,2},trans,dynamics{1,2});
loc_ = loc;
loc_(2) = location('loc2',inv{1,2},trans,dynamics_);

% init first hybrid automaton
HA1 = hybridAutomaton(loc);
HA1_ = hybridAutomaton(loc_);

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
reset_ = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[-1;-1;1]);
trans = transition(guard,reset,1);
trans_ = transition(guard,reset_,1); 
loc(2) = location('loc2',inv{2,2},trans,dynamics{2,2});
loc_ = loc;
loc_(2) = location('loc2',inv{2,2},trans_,dynamics{2,2});

% init second hybrid automaton
HA2 = hybridAutomaton(loc);
HA2_ = hybridAutomaton(loc_);

% components and input binds
components = [HA1;HA2];
components_ = [HA1_;HA2];
components__ = [HA1;HA2_];
inputBinds{1} = [2 1; 2 2];
inputBinds{2} = [1 1];

% instantiate parallel hybrid automaton
pHA = parallelHybridAutomaton(components,inputBinds);
% different flow
pHA_ = parallelHybridAutomaton(components_,inputBinds);
% different reset function
pHA__ = parallelHybridAutomaton(components__,inputBinds);

% full comparison
res(end+1,1) = isequal(pHA,pHA);
temp = isequal([pHA,pHA_],[pHA,pHA__]);
res(end+1,1) = all(size(temp) == [1,2]);
res(end+1,1) = all(temp == [true false]);
% slightly different automata
res(end+1,1) = ~isequal(pHA,pHA_);
res(end+1,1) = ~isequal(pHA,pHA__);
res(end+1,1) = ~isequal([pHA,pHA_],pHA);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
