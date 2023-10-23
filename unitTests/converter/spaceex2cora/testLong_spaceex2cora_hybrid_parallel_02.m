function res = testLong_spaceex2cora_hybrid_parallel_02
% testLong_spaceex2cora_hybrid_parallel_02 - test for model
%    conversion from SpaceEx to CORA for a simple parallel hybrid automaton
%    with two components with all-nonlinear behavior
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_parallel_02
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       12-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;


% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% file name of SpaceEx model file
filename = 'test_hybrid_parallel_twocomp2';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model
syms x y z

% first component, first location
f = @(x,u) [-sin(x(2)) + u(1); x(1)^2 + u(2)];
g = @(x,u) 0.05*(x(1)^2 + x(2));
dynamics = nonlinearSys([filename '_Comp1_Loc1_FlowEq'],f,g);
eq = x + y^2;
inv = levelSet(eq,[x;y],'<=');
guard = levelSet([eq,-eq],[x;y],{'<=','<='});
reset.f = @(x,u) [sqrt(x(1)) + 3; x(2)^2 + 3];
trans = transition(guard,reset,2);
loc(1) = location('loc1',inv,trans,dynamics);

% first component, second location
f = @(x,u) [-x(2)^3 + u(2); -cos(x(1)) + u(1)];
g = @(x,u) 0.05*(x(1) - x(2)^2);
dynamics = nonlinearSys([filename '_Comp1_Loc2_FlowEq'],f,g);
eq = -x - y^2;
inv = levelSet(eq,[x;y],'<=');
guard = levelSet([eq,-eq],[x;y],{'<=','<='});
reset.f = @(x,u) [x(1)^2 - 3; x(2)^3 + 3];
trans = transition(guard,reset,1);
loc(2) = location('loc2',inv,trans,dynamics);

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
f = @(x,u) [x(2) - x(3)^2; cos(x(1)) - u(1); x(2)^2];
g = @(x,u) [0.05*sin(x(3)); 0.05*(x(1) + cos(x(2)))];
dynamics = nonlinearSys([filename '_Comp2_Loc1_FlowEq'],f,g);
eq = x^2 + y + z - 1;
inv = levelSet(eq,[x;y;z],'<=');
guard = levelSet([eq,-eq],[x;y;z],{'<=','<='});
reset.f = @(x,u) [x(1)^2 + 1; x(2) + 1; x(3) + 1];
trans = transition(guard,reset,2);
loc(1) = location('loc1',inv,trans,dynamics);

% second component, second location
f = @(x,u) [-x(2)^2 + x(3); sin(x(1)); cos(x(2)) + u(1)];
g = @(x,u) [0.05*sin(x(1)); 0.05*(x(2) - cos(x(3)))];
dynamics = nonlinearSys([filename '_Comp2_Loc2_FlowEq'],f,g);
eq = -x^2 - y - z + 1;
inv = levelSet(eq,[x;y;z],'<=');
guard = levelSet([eq,-eq],[x;y;z],{'<=','<='});
reset.f = @(x,u) [x(1)^2 - 1; x(2) - 1; x(3) - 1];
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
