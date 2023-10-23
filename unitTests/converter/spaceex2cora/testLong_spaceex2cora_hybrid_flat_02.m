function res = testLong_spaceex2cora_hybrid_flat_02
% testLong_spaceex2cora_hybrid_flat_02 - test for model conversion 
%    from SpaceEx to CORA for a simple hybrid system with two locations
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_flat_02
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
filename = 'test_hybrid_flat_twoloc1';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model

% loc1
inv = polytope([-1 0],0);

% transitions
c = [1;0]; d = 0;
guard = conHyperplane(c,d);
reset = struct('A',[0,0;0,1],'c',[5;0]);
trans = transition(guard,reset,2);

% flow equation
f = @(x,u) [-2*sin(x(1)) + u(1); x(1) - x(2)];
dynamics = nonlinearSys([filename '_Loc1_FlowEq'],f);

% define location
loc(1) = location('loc1',inv,trans,dynamics);

% loc2
inv = polytope([1 0],0);

% transitions
c = [1;0]; d = 0;
guard = conHyperplane(c,d);
reset = struct('A',[0,0;0,1],'c',[-5;0]);
trans = transition(guard,reset,1);

% flow equation
f = @(x,u) [-2*cos(x(1)) + u(1); x(1) + x(2)];
dynamics = nonlinearSys([filename '_Loc2_FlowEq'],f);

% define location
loc(2) = location('loc2',inv,trans,dynamics);

% instantiate hybrid automaton
sys_cora = hybridAutomaton(loc);


% compare systems
if sys_cora ~= sys_spaceex
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
