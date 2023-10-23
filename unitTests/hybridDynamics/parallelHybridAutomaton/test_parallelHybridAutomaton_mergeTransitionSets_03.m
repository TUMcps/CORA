function res = test_parallelHybridAutomaton_mergeTransitionSets_03
% test_parallelHybridAutomaton_mergeTransitionSets_03 - test function for
%    merging transitions in the location product
%
% Syntax:
%    res = test_parallelHybridAutomaton_mergeTransitionSets_03
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
% Written:       24-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% first component, first location
dynamics{1,1} = linearSys([0 -1; 1 0],[1 0; 0 1],[],[0.05 0.05]);
inv{1,1} = polytope([1 1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[3;3]);
trans{1,1} = transition(guard,reset,2,'sync1');
loc = location('loc1',inv{1,1},trans{1,1},dynamics{1,1});

% first component, second location
dynamics{1,2} = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
inv{1,2} = polytope([-1 -1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[-3;3]);
trans{1,2} = transition(guard,reset,1,'sync2');
loc(2) = location('loc2',inv{1,2},trans{1,2},dynamics{1,2});

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
dynamics{2,1} = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
inv{2,1} = polytope([1 1 1],1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[1;1;1]);
trans{2,1} = transition(guard,reset,2,'sync1');
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

% third component, first location
dynamics{3,1} = linearSys([0 1; 1 0],[0;1],[],[0 0.05; 0.05 0]);
inv{3,1} = polytope([1 1],1);
guard = conHyperplane([1 1],1);
reset = struct('A',[1 0; 0 1],'c',[1;1]);
trans{3,1} = transition(guard,reset,2,'sync2');
loc(1) = location('loc1',inv{3,1},trans{3,1},dynamics{3,1});

% third component, second location
dynamics{3,2} = linearSys([0 -1; 1 0],[0;1],[],[0.05 0; 0 -0.05]);
inv{3,2} = polytope([-1 -1],-1);
guard = conHyperplane([1 1],1);
reset = struct('A',[1 0; 0 1],'c',[-1;-1]);
trans{3,2} = transition(guard,reset,1);
loc(2) = location('loc2',inv{3,2},trans{3,2},dynamics{3,2});

% init second hybrid automaton
HA3 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2;HA3];
inputBinds{1} = [2 1; 2 2];
inputBinds{2} = [1 1];
inputBinds{3} = [0 1];

% instantiate parallel hybrid automaton
pHA = parallelHybridAutomaton(components,inputBinds);

% all synchronization labels incl. information
allLabels = struct('name',{'sync1','sync2'},...
    'component',{[1;2],[1;3]},...
    'location',{[1;1],[2;1]},...
    'transition',{[1;1],[1;1]},...
    'guardempty',{[false;false],[false;false]});

% merge invariants of locations [1;1;2]
% comp1,loc1,trans1:
% - guard = conHyperplane([1 1],0);
% - reset = struct('A',[1 0; 0 1],'c',[3;3]);
% - trans{1,1}{1} = transition(guard,reset,2,'sync1');
% comp2,loc1,trans1:
% - guard = conHyperplane([1 1 1],1);
% - reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[1;1;1]);
% - trans{2,1}{1} = transition(guard,reset,2,'sync1');
% comp3,loc2,trans1:
% - guard = conHyperplane([1 1],1);
% - reset = struct('A',[1 0; 0 1],'c',[-1;-1]);
% - trans{3,2}{1} = transition(guard,reset,1);

locID = [1;1;2];

% merge transition sets
mergedTrans = mergeTransitionSets(pHA,...
    {pHA.components(1).location(locID(1)).transition,...
    pHA.components(2).location(locID(2)).transition,...
    pHA.components(3).location(locID(3)).transition},...
    locID,allLabels);

% manual merge of transitions:
% transitions of comp1/loc1 and comp2/loc1 have the same synchronization
% label, thus can only occur together; additionally, comp3/loc2 has no
% synchronization label and can thus occur whenever -> two transitions

% first transition: 'sync1' label
target = [2;2;2];

% guard set
guard = lift(polytope(trans{1,locID(1)}.guard),7,1:2) ...
    & lift(polytope(trans{2,locID(2)}.guard),7,3:5);

% reset function: identity for third component
reset.A = blkdiag(trans{1,locID(1)}.reset.A,...
    trans{2,locID(2)}.reset.A,...
    eye(length(trans{3,locID(3)}.reset.A)));
reset.c = [trans{1,locID(1)}.reset.c;...
    trans{2,locID(2)}.reset.c;...
    zeros(length(trans{3,locID(3)}.reset.c),1)];

% full first transition
mergedTrans_(1) = transition(guard,reset,target);


% second transition: only third component
target = [1;1;1];

% guard set
guard = lift(trans{3,locID(3)}.guard,7,6:7);

% reset function: identity for first/second component
reset.A = blkdiag(eye(length(trans{1,locID(1)}.reset.A)),...
    eye(length(trans{2,locID(2)}.reset.A)),...
    trans{3,locID(3)}.reset.A);
reset.c = [zeros(length(trans{1,locID(1)}.reset.c),1);...
    zeros(length(trans{2,locID(2)}.reset.c),1);...
    trans{3,locID(3)}.reset.c];

% full second transition
mergedTrans_(2) = transition(guard,reset,target);


% compare solutions
res(end+1,1) = length(mergedTrans) == length(mergedTrans_);
res(end+1,1) = isequal(mergedTrans(1),mergedTrans_(1),1e-14);
res(end+1,1) = isequal(mergedTrans(2),mergedTrans_(2),1e-14);


% merge invariants of locations [2;1;2]
% comp1,loc2,trans1:
% - guard = conHyperplane([1 1],0);
% - reset = struct('A',[1 0; 0 1],'c',[-3;3]);
% - trans{1,2}{1} = transition(guard,reset,1,'sync2');
% comp2,loc1,trans1:
% - guard = conHyperplane([1 1 1],1);
% - reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[1;1;1]);
% - trans{2,1}{1} = transition(guard,reset,2,'sync1');
% comp3,loc2,trans1:
% - guard = conHyperplane([1 1],1);
% - reset = struct('A',[1 0; 0 1],'c',[-1;-1]);
% - trans{3,2}{1} = transition(guard,reset,1);

locID = [2;1;2];

% merge transition sets
mergedTrans = mergeTransitionSets(pHA,...
    {pHA.components(1).location(locID(1)).transition,...
    pHA.components(2).location(locID(2)).transition,...
    pHA.components(3).location(locID(3)).transition},...
    locID,allLabels);

% manual merge of transitions:
% sync labels 'sync1' and 'sync2' cannot trigger, because other components
% are not in correct location (with these labels); hence, only the
% transition from component 3 can be taken

% transition from loc2 to loc1 (of component 3)
target = [2;1;1];

% guard set
guard = lift(trans{3,locID(3)}.guard,7,6:7);

% reset function: identity for third component
reset.A = blkdiag(eye(length(trans{1,locID(1)}.reset.A)),...
    eye(length(trans{2,locID(2)}.reset.A)),...
    trans{3,locID(3)}.reset.A);
reset.c = [zeros(length(trans{1,locID(1)}.reset.c),1);...
    zeros(length(trans{2,locID(2)}.reset.c),1);...
    trans{3,locID(3)}.reset.c];

% full first transition
mergedTrans_ = transition(guard,reset,target);


% compare solutions
res(end+1,1) = length(mergedTrans) == length(mergedTrans_);
res(end+1,1) = isequal(mergedTrans,mergedTrans_,1e-14);

% combine results
res = all(res);


% ------------------------------ END OF CODE ------------------------------
