function res = test_parallelHybridAutomaton_mergeTransitionSets_01
% test_parallelHybridAutomaton_mergeTransitionSets_01 - test function for
%    merging transitions in the location product
%
% Syntax:
%    res = test_parallelHybridAutomaton_mergeTransitionSets_01
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
% Written:       15-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% first component, first location
dynamics{1,1} = linearSys([0 -1; 1 0],[1 0; 0 1],[],[0.05 0.05]);
inv{1,1} = polytope([1 1],0);
guard = polytope([],[],[1 1],0);
reset = linearReset([1 0; 0 1],zeros(2,2),[3;3]);
trans{1,1} = transition(guard,reset,2);
loc = location('loc1',inv{1,1},trans{1,1},dynamics{1,1});

% first component, second location
dynamics{1,2} = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
inv{1,2} = polytope([-1 -1],0);
guard = polytope([],[],[1 1],0);
reset = linearReset([1 0; 0 1],zeros(2,2),[-3;3]);
trans{1,2} = transition(guard,reset,1);
loc(2) = location('loc2',inv{1,2},trans{1,2},dynamics{1,2});

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
dynamics{2,1} = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
inv{2,1} = polytope([1 1 1],1);
guard = polytope([],[],[1 1 1],1);
reset = linearReset([1 0 0; 0 1 0; 0 0 1],zeros(3,1),[1;1;1]);
trans{2,1} = transition(guard,reset,2);
loc(1) = location('loc1',inv{2,1},trans{2,1},dynamics{2,1});

% second component, second location
dynamics{2,2} = linearSys([0 -1 1; 1 0 0; 0 1 0],[0;0;1],[],[0.05 0 0; 0 0.05 -0.05]);
inv{2,2} = polytope([-1 -1 -1],-1);
guard = polytope([],[],[1 1 1],1);
reset = linearReset([1 0 0; 0 1 0; 0 0 1],zeros(3,1),[-1;-1;-1]);
trans{2,2} = transition(guard,reset,1);
loc(2) = location('loc2',inv{2,2},trans{2,2},dynamics{2,2});

% init second hybrid automaton
HA2 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2];
inputBinds{1} = [2 1; 2 2];
inputBinds{2} = [1 1];

% instantiate parallel hybrid automaton
pHA = parallelHybridAutomaton(components,inputBinds);

% no labels in transitions
allLabels = struct([]);

% merge invariants of different location pairs (first column: location of
% component 1, second column: location of component 2)
comb = [1 1; 1 2; 2 1; 2 2];

for i=1:size(comb,1)

    % merge transition sets
    mergedTrans = mergeTransitionSets(pHA,comb(i,:)',allLabels);
    
    % instantiate true solution: 2 transitions

    % guard set
    guard = polytope([],[],[trans{1,comb(i,1)}.guard.Ae,...
        zeros(1,length(trans{2,comb(i,2)}.guard.Ae))],0);

    % reset function
    reset = linearReset(eye(5),[],...
        [trans{1,comb(i,1)}.reset.c;...
        zeros(length(trans{2,comb(i,2)}.reset.c),1)]);

    % new target location ID
    target = [0;comb(i,2)];
    if comb(i,1) == 1; target(1) = 2; else; target(1) = 1; end

    % full first transition
    mergedTrans_(1) = transition(guard,reset,target);

    % guard set
    guard = polytope([],[],[zeros(1,length(trans{1,comb(i,1)}.guard.Ae)),...
        trans{2,comb(i,2)}.guard.Ae],1);

    % reset function
    reset = linearReset(eye(5),[],...
        [zeros(length(trans{1,comb(i,1)}.reset.c),1);...
        trans{2,comb(i,2)}.reset.c]);
    
    % new target location ID
    target = [comb(i,1);0];
    if comb(i,2) == 1; target(2) = 2; else; target(2) = 1; end

    % full second transition
    mergedTrans_(2) = transition(guard,reset,target);

    % compare solutions
    assert(length(mergedTrans) == length(mergedTrans_));
    assert(isequal(mergedTrans(1),mergedTrans_(1),1e-14));
    assert(isequal(mergedTrans(2),mergedTrans_(2),1e-14));
end

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
