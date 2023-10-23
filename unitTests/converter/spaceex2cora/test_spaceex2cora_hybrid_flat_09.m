function res = test_spaceex2cora_hybrid_flat_09
% test_spaceex2cora_hybrid_flat_09 - test for model conversion from
%    SpaceEx to CORA for a simple hybrid automaton with two
%    components with empty invariants but given guard sets; the conversion
%    then has to compute the invariants to match the respective guard sets
%
% Syntax:
%    test_spaceex2cora_hybrid_flat_09
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       25-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;


% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% file name of SpaceEx model file
filename = 'test_hybrid_flat_twoloc4';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model

% loc1
inv = emptySet(2);

% transitions
syms x y
eq = x^2 + y^2 - 5;
guard = levelSet(eq,[x;y],'<=');
reset = struct('A',[1,0;0,1],'c',[-2;-2]);
trans(1) = transition(guard,reset,2);
eq = -x^2 - y^2 + 5;
guard = levelSet(eq,[x;y],'<=');
reset = struct('A',[1,0;0,1],'c',[2;2]);
trans(2) = transition(guard,reset,2);

% flow equation
f = @(x,u) [-2*sin(x(1)) + u(1); x(1) - x(2)];
dynamics = nonlinearSys([filename '_Loc1_FlowEq'],f);

% define location
loc(1) = location('loc1',inv,trans,dynamics);

% loc2
inv = fullspace(2);

% no transition
trans = transition();

% flow equation
f = @(x,u) [-2*cos(x(1)) + u(1); x(1) + x(2)];
dynamics = nonlinearSys([filename '_Loc2_FlowEq'],f);

% define location
loc(2) = location('loc2',inv,trans,dynamics);

% instantiate hybrid automaton
sys_cora = hybridAutomaton(loc);


% compare systems
if ~isequal(sys_cora,sys_spaceex,1e-4)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
