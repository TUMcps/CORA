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

% assume true
res = [];

% first component, first location
dynamics{1,1} = linearSys([0 -1; 1 0],[1 0; 0 1],[],[0.05 0.05]);
inv{1,1} = polytope([1 1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[3;3]);
trans{1,1} = transition(guard,reset,2);
loc = location('loc1',inv{1,1},trans{1,1},dynamics{1,1});

% first component, second location
dynamics{1,2} = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
inv{1,2} = polytope([-1 -1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[-3;3]);
trans{1,2} = transition(guard,reset,1);
loc(2) = location('loc2',inv{1,2},trans{1,2},dynamics{1,2});

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
dynamics{2,1} = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
inv{2,1} = polytope([1 1 1],1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[1;1;1]);
trans{2,1} = transition(guard,reset,2);
loc(1) = location('loc1',inv{2,1},trans{2,1},dynamics{2,1});

% second component, second location
dynamics{2,2} = linearSys([0 -1 1; 1 0 0; 0 1 0],[0;0;1],[],[0.05 0 0; 0 0.05 -0.05]);
inv{2,2} = polytope([-1 -1 -1],-1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[-1;-1;-1]);
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
    mergedTrans = mergeTransitionSets(pHA,...
        {pHA.components(1).location(comb(i,1)).transition,...
        pHA.components(2).location(comb(i,2)).transition},...
        comb(i,:)',allLabels);
    
    % instantiate true solution: 2 transitions

    % guard set
    guard = conHyperplane([trans{1,comb(i,1)}.guard.a,...
        zeros(1,length(trans{2,comb(i,2)}.guard.a))],0);

    % reset function
    reset = struct('A',eye(5),'c',...
        [trans{1,comb(i,1)}.reset.c;...
        zeros(length(trans{2,comb(i,2)}.reset.c),1)]);

    % new target location ID
    target = [0;comb(i,2)];
    if comb(i,1) == 1; target(1) = 2; else; target(1) = 1; end

    % full first transition
    mergedTrans_(1) = transition(guard,reset,target);

    % guard set
    guard = conHyperplane([zeros(1,length(trans{1,comb(i,1)}.guard.a)),...
        trans{2,comb(i,2)}.guard.a],1);

    % reset function
    reset = struct('A',eye(5),'c',...
        [zeros(length(trans{1,comb(i,1)}.reset.c),1);...
        trans{2,comb(i,2)}.reset.c]);
    
    % new target location ID
    target = [comb(i,1);0];
    if comb(i,2) == 1; target(2) = 2; else; target(2) = 1; end

    % full second transition
    mergedTrans_(2) = transition(guard,reset,target);

    % compare solutions
    res(end+1,1) = length(mergedTrans) == length(mergedTrans_);
    res(end+1,1) = isequal(mergedTrans(1),mergedTrans_(1),1e-14);
    res(end+1,1) = isequal(mergedTrans(2),mergedTrans_(2),1e-14);
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
