function res = test_spaceex2cora_hybrid_flat_04
% test_spaceex2cora_hybrid_flat_04 - test for model conversion from SpaceEx
%    to CORA for a simple hybrid system with four locations
%
% Syntax:
%    test_spaceex2cora_hybrid_flat_04
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
filename = 'test_hybrid_flat_fourloc';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model

% top-left
inv = polytope([1 0; 0 -1],[0; 0]);
dynamics = linearSys([0 0; 0 0],[1; 0],[0; -1]);
% transitions
guard = conHyperplane([-1 0],0,[0 -1],0);
reset = struct('A',[1,0;0,1],'c',[10;5]);
trans(1) = transition(guard,reset,2);
guard = conHyperplane([0 1],0,[-1 0],0);
reset = struct('A',[1,0;0,1],'c',[-5;-10]);
trans(2) = transition(guard,reset,3);
% define location
loc = location('topleft',inv,trans,dynamics);

% top-right
inv = polytope([-1 0; 0 -1],[0; 0]);
dynamics = linearSys([0 0; 0 0],[0; 1],[-1; 0]);
% transitions
guard = conHyperplane([1 0],0,[0 -1],0);
reset = struct('A',[1,0;0,1],'c',[-10;5]);
trans(1) = transition(guard,reset,1);
guard = conHyperplane([0 1],0,[-1 0],0);
reset = struct('A',[1,0;0,1],'c',[5;-10]);
trans(2) = transition(guard,reset,4);
% define location
loc(2) = location('topright',inv,trans,dynamics);

% bottom-left
inv = polytope([1 0; 0 1],[0; 0]);
dynamics = linearSys([0 0; 0 0],[0; 1],[1; 0]);
% transitions
guard = conHyperplane([0 -1],0,[1 0],0);
reset = struct('A',[1,0;0,1],'c',[-5;10]);
trans(1) = transition(guard,reset,1);
guard = conHyperplane([-1 0],0,[0 1],0);
reset = struct('A',[1,0;0,1],'c',[10;-5]);
trans(2) = transition(guard,reset,4);
% define location
loc(3) = location('bottomleft',inv,trans,dynamics);

% bottom-right
inv = polytope([-1 0; 0 1],[0; 0]);
dynamics = linearSys([0 0; 0 0],[1; 0],[0; 1]);
% transitions
guard = conHyperplane([0 -1],0,[-1 0],0);
reset = struct('A',[1,0;0,1],'c',[5;10]);
trans(1) = transition(guard,reset,2);
guard = conHyperplane([1 0],0,[0 1],0);
reset = struct('A',[1,0;0,1],'c',[-10;-5]);
trans(2) = transition(guard,reset,3);
% define location
loc(4) = location('bottomright',inv,trans,dynamics);

% instantiate hybrid automaton
sys_cora = hybridAutomaton(loc);


% compare systems
if sys_cora ~= sys_spaceex
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
