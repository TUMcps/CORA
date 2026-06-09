function res = test_cora2spaceex_hybrid
% test_cora2spaceex_hybrid - test for model conversion from CORA to SpaceEx
%    for a hybrid system
%
% Syntax:
%    test_cora2spaceex_hybrid
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Niklas Kochdumper
% Written:       14-November-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% load model
HA = roomHeating();

% rotate state space of the hybrid automaton (to generate non-integer
% double numbers in all components)
phi = rand()*2*pi; c = rand(2,1);
R = [cos(phi) -sin(phi); sin(phi) cos(phi)];

locs = [];

for i = 1:length(HA.location)

    loc = HA.location(i);

    % modify system dynamcis
    inva = R*loc.invariant + c;
    sys = loc.contDynamics;
    sys = linearSys(R*sys.A,R*sys.B,R*sys.c);

    trans = [];

    % modify all transitions
    for j = 1:length(loc.transition)
        guard = (R*loc.transition(j).guard + c) & inva;
        reset = loc.transition(j).reset;
        reset = linearReset(R*reset.A,[],R*reset.c);
        trans = [trans;transition(guard,reset,loc.transition(j).target)];
    end

    locs = [locs;location(inva,trans,sys)];
end

HA = hybridAutomaton(locs);

% convert model to SpaceEx format
cora2spaceex(HA,'model_test_cora2spaceex_hybrid');

% import model from SpaceEx format
spaceex2cora('model_test_cora2spaceex_hybrid',[],[],...
                                    'conv_test_cora2spaceex_hybrid');
HA_ = conv_test_cora2spaceex_hybrid();

% compare hybrid automata
assert(isequal(HA,HA_,100*eps));

% ------------------------------ END OF CODE ------------------------------
