function res = test_spaceex2cora_hybrid_flat_03
% test_spaceex2cora_hybrid_flat_03 - test for model conversion from SpaceEx
%    to CORA for a simple hybrid system with one location
%
% Syntax:
%    test_spaceex2cora_hybrid_flat_03
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

% assume true
res = true;


% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% file name of SpaceEx model file
filename = 'test_hybrid_flat_oneloc3';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model
inv = polytope([-1 0; 0 1],[0; 0]);

% transitions
c = [1;0]; d = 0; C = [0 1]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',[-1,0;0,1],'B',[0;1],'c',[0;0]);
trans(1) = transition(guard,reset,1);
c = [0;-1]; d = 0; C = [1 0]; D = 0;
guard = conHyperplane(c,d,C,D);
reset = struct('A',[0,1;-1,0],'c',[-1;0]);
trans(2) = transition(guard,reset,1);

% flow equation
dynamics = linearSys([-2 1;1 -2],[-1;0]);

% define location
loc = location('always',inv,trans,dynamics);

% instantiate hybrid automaton
sys_cora = hybridAutomaton(loc);

% compare systems
if sys_cora ~= sys_spaceex
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
