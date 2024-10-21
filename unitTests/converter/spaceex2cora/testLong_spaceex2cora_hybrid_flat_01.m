function res = testLong_spaceex2cora_hybrid_flat_01
% testLong_spaceex2cora_hybrid_flat_01 - test for model conversion
%    from SpaceEx to CORA for a simple hybrid system with one location
%
% Syntax:
%    res = testLong_spaceex2cora_hybrid_flat_01
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       10-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% file name of SpaceEx model file
filename = 'test_hybrid_flat_oneloc4';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model
inv = polytope([-1 0; 0 1],[0; 0]);

% transitions
guard = polytope([0,1],0,[1,0],0);
reset = nonlinearReset(@(x,u) [-x(1); sin(x(2)) + u(1)],2,1,2);
trans(1) = transition(guard,reset,1);

guard = polytope([-1,0],0,[0,-1],0);
reset = linearReset([0,1;-1,0],[0;0],[-1;0]);
trans(2) = transition(guard,reset,1);

% flow equation
f = @(x,u) [-2*sin(x(1)) + u(1); log(x(1)) - x(2)];
dynamics = nonlinearSys([filename '_Loc1_FlowEq'],f);

% define location
loc = location('always',inv,trans,dynamics);

% instantiate hybrid automaton
sys_cora = hybridAutomaton(loc);

% compare systems
assert(sys_cora == sys_spaceex);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
