function res = testLong_spaceex2cora_hybrid_parallel_01
% testLong_spaceex2cora_hybrid_parallel_01 - test for model
%    conversion from SpaceEx to CORA for a simple parallel hybrid automaton
%    with two components with all-linear behavior
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_parallel_01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       11-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;


% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% file name of SpaceEx model file
filename = 'test_hybrid_parallel_twocomp1';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model

% first component, first location
dynamics = linearSys([0 -1; 1 0],[1 0; 0 1],[],[0.05 0.05]);
inv = polytope([1 1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[3;3]);
trans = transition(guard,reset,2);
loc(1) = location('loc1',inv,trans,dynamics);

% first component, second location
dynamics = linearSys([0 -1; -1 0],[0 1; 1 0],[],[0.05 -0.05]);
inv = polytope([-1 -1],0);
guard = conHyperplane([1 1],0);
reset = struct('A',[1 0; 0 1],'c',[-3;3]);
trans = transition(guard,reset,1);
loc(2) = location('loc2',inv,trans,dynamics);

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
dynamics = linearSys([0 1 -1; 1 0 0; 0 1 0],[0;-1;0],[],[0 0 0.05; 0.05 0.05 0]);
inv = polytope([1 1 1],1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[1;1;1]);
trans = transition(guard,reset,2);
loc(1) = location('loc1',inv,trans,dynamics);

% second component, second location
dynamics = linearSys([0 -1 1; 1 0 0; 0 1 0],[0;0;1],[],[0.05 0 0; 0 0.05 -0.05]);
inv = polytope([-1 -1 -1],-1);
guard = conHyperplane([1 1 1],1);
reset = struct('A',[1 0 0; 0 1 0; 0 0 1],'c',[-1;-1;-1]);
trans = transition(guard,reset,1);
loc(2) = location('loc2',inv,trans,dynamics);

% init second hybrid automaton
HA2 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2];
inputBinds{1} = [2 1; 2 2];
inputBinds{2} = [1 1];

% instantiate parallel hybrid automaton
sys_cora = parallelHybridAutomaton(components,inputBinds);


% compare systems
if ~isequal(sys_cora,sys_spaceex,1e-4)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
