function res = test_linearSys_homogeneousSolution
% test_linearSys_homogeneousSolution - unit test for the computation of the
%    homogeneous time-point and time-interval solution
%
% Syntax:
%    res = test_linearSys_homogeneousSolution
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
% Written:       16-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-14;

% init system, state set, and algorithm parameters
A = [-1 -4; 4 -1];
sys = linearSys(A);
X = zonotope([40;20],[1 4 2; -1 3 5]);
timeStep = 0.01;
truncationOrder = 5;

% compute homogeneous solutions and error term
[Htp,Hti,C_state] = homogeneousSolution(sys,X,timeStep,truncationOrder);

% time-point solution: e^Adt * X
eAdt = expm(A*timeStep);
assert(isequal(Htp,eAdt*X,tol));

% time-interval solution must contain all e^At*x with t in [0,dt], x in X;
% we choose the vertices for this
V = vertices(X);
t = linspace(0,timeStep,101);
eAt_V = cell2mat(arrayfun(@(t_) expm(A*t_)*V,t,'UniformOutput',false));
assert(all(contains(Hti,eAt_V,'exact',tol)));

% the error term is computed with the following meaning: the true state
% e^At*x with t in [0,dt], x in X must be contained in the linear
% interpolation between x in X and e^At*x at t=dt enlarged by the error
% term
for i=1:100
    % sample state and time
    t = rand*timeStep;
    x = randPoint(X);

    % true propagated state
    x_prop_true = expm(A*t)*x;

    % linear interpolation
    eAdt_x = eAdt*x;
    x_prop_lint = x + t/timeStep * (eAdt_x - x);

    % check containment   
    assertLoop(contains(x_prop_lint + C_state,x_prop_true,'exact',tol), i);
end


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
