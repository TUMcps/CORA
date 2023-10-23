function res = test_parallelHybridAutomaton_mergeInvariants_01
% test_parallelHybridAutomaton_mergeInvariants_01 - test function for
%    merging invariant sets (polytopes) in the location product
%
% Syntax:
%    res = test_parallelHybridAutomaton_mergeInvariants_01
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
res = true;

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


% merge invariants of different location pairs
comb = [1 1; 1 2; 2 1; 2 2];

for i=1:size(comb,1)

    % merge invariants
    mergedInv = mergeInvariants(pHA,...
        {pHA.components(1).location(comb(i,1)).invariant,...
        pHA.components(2).location(comb(i,2)).invariant},...
        {pHA.components(1).location(comb(i,1)).transition,...
        pHA.components(2).location(comb(i,2)).transition},{});
    
    % instantiate true solution
    mergedInv_ = polytope(blkdiag(inv{1,comb(i,1)}.A,inv{2,comb(i,2)}.A),...
        [inv{1,comb(i,1)}.b;inv{2,comb(i,2)}.b]);
    
    % compare solutions
    if ~eq(mergedInv,mergedInv_,1e-14)
        res = false;
    end
end

% ------------------------------ END OF CODE ------------------------------
