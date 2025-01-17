function res = test_linearSys_reachBackward_pursuitevasion
% test_linearSys_reachBackward_pursuitevasion - unit test for backward
%    reachability analysis of a simple benchmark, checks approximate
%    containment relations for a one-step time-interval backward reachable
%    set
%
% Syntax:
%    res = test_linearSys_reachBackward_pursuitevasion
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       20-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% system dynamics: relative position and velocity of two quadrotors
A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0];
% input (player 1): relative acceleration in x and y
B = [0 0;
     1 0;
     0 0;
     0 1];
% disturbance (player 2): relative acceleration in x and y
E = [0  0;
    -1  0;
     0  0;
     0 -1];
sys = linearSys('sys',A,B,[],[],[],[],E);

% times
params.tStart = 0;
params.tFinal = 0.05;
options.timeStep = 0.05;

% target set = set of collisions: relative position between quadrotors
params.R0 = polytope(interval(-ones(4,1),ones(4,1)));
% set of controllable inputs: acceleration in x and y of player 1
U = zonotope(zeros(2,1),eye(2));
% set of uncontrollable disturbances: acceleration in x and y of player 2
W = zonotope(zeros(2,1),eye(2));

% EA backward reachability
options.linAlg = 'inner:EA:timeinterval';
% - no inputs, with disturbances
params.U = 0*U;
params.W = W;
R_small = reachBackward(sys,params,options);
% - no inputs, no disturbances
params.U = 0*U;
params.W = 0*W;
R_medium = reachBackward(sys,params,options);
% - with inputs, no disturbances
params.U = U;
params.W = 0*W;
R_large = reachBackward(sys,params,options);

% the more input capacity vs. disturbance capacity, the larger the backward
% reachable set must be
assert(aux_checkContainment(R_large.timeInterval.set{end},...
                            R_medium.timeInterval.set{end}));
assert(aux_checkContainment(R_medium.timeInterval.set{end},...
                            R_small.timeInterval.set{end}));


% AE backward reachability
options.linAlg = 'outer:AE:timeinterval';
% - with inputs, no disturbances
params.U = U;
params.W = 0*W;
R_small = reachBackward(sys,params,options);
% - no inputs, no disturbances
params.U = 0*U;
params.W = 0*W;
R_medium = reachBackward(sys,params,options);
% - no inputs, with disturbances
params.U = 0*U;
params.W = W;
R_large = reachBackward(sys,params,options);

% the more disturbance capacity vs. input capacity, the larger the backward
% reachable set must be
assert(aux_checkContainment(R_large.timeInterval.set{end},...
                            R_medium.timeInterval.set{end}));
assert(aux_checkContainment(R_medium.timeInterval.set{end},...
                            R_small.timeInterval.set{end}));


% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkContainment(S1,S2)
% use support function evaluations in axis-aligned and diagonal directions
% to check whether S1 contains S2 (both are convex)

tol = 1e-8;
n = dim(S1);
directions = [eye(n), -eye(n), vertices(interval(-ones(n,1),ones(n,1)))];

for i=1:size(directions,2)
    val1 = supportFunc_(S1,directions(:,i),'upper');
    val2 = supportFunc_(S2,directions(:,i),'upper');
    if val1 < val2 && ~withinTol(val1,val2,tol)
        res = false;
        return
    end
end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
