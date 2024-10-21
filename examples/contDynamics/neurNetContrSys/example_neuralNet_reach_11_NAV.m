function [completed,res,tTotal] = example_neuralNet_reach_11_NAV(type)
% example_neuralNet_reach_11_NAV - example of reachability analysis
%    for an neural network controlled system
%
% Syntax:
%    completed = example_neuralNet_reach_11_NAV()
%    completed = example_neuralNet_reach_11_NAV(type)
%
% Inputs:
%    type - str, 'set' and 'point'
%
% Outputs:
%    completed - boolean
%    res - verification result
%    tTotal - total time
%

% Authors:       Manuel Wendl, Tobias Ladner
% Written:       30-April-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin == 0
    type = 'set';
end
inputArgsCheck({{type,'str',{'point','set'}}})

fprintf("BENCHMARK: NavigationTask (NAV-%s)\n", type)

% Parameter ---------------------------------------------------------------

params.tFinal = 6;
params.R0 = zonotope(interval( ...
    [2.9; 2.9; 0; 0], ...
    [3.1; 3.1; 0; 0] ...
));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 1;
options.zonotopeOrder = 80;

% Parameters for NN evaluation --------------------------------------------

options.nn.poly_method = 'regression';

% System Dynamics ---------------------------------------------------------

% open-loop system
f = @(x,u) [x(3)*cos(x(4));x(3)*sin(x(4));u(1);u(2)];
sys = nonlinearSys(f);

% load neural network controller
if strcmp(type,'point')
    nn = neuralNetwork.readONNXNetwork('nn-nav-point.onnx');
else
    nn = neuralNetwork.readONNXNetwork('nn-nav-set.onnx');
end

% construct neural network controlled system
sys = neurNetContrSys(sys, nn, 0.2);

% Specification -----------------------------------------------------------

specGoalSet = specification(interval( ...
    [-0.5;-0.5;-inf;-inf], ...
    [0.5;0.5;inf;inf] ...
), 'safeSet', interval(params.tFinal));

specObstacle = specification(interval( ...
    [1;1;-inf;-inf], ...
    [2;2;inf;inf] ...
), 'unsafeSet');

spec = add(specGoalSet,specObstacle);

% Verification ------------------------------------------------------------

t = tic;
[res, R, simRes] = verify(sys, spec, params, options, true);
tTotal = toc(t);
disp(['Result: ' res])

% Visualization -----------------------------------------------------------
disp("Plotting..")

figure; hold on; box on;

% plot specification (over entire time horizon)
plot(specGoalSet, [1,2], 'DisplayName', 'Goal set');
plot(specObstacle, [1,2], 'DisplayName', 'Obstacle');

% plot reachable set
useCORAcolors('CORA:contDynamics')
plot(R, [1,2], 'Unify', true, 'DisplayName', 'Reachable set');

% plot initial set
plot(R(1).R0,[1,2],'DisplayName','Initial set');

% plot simulation
plot(simRes, [1,2], 'DisplayName', 'Simulations','Color','k');

% labels and legend
xlabel('x'); ylabel('y');
xlim([-1, 3.5]); ylim([-1, 4]);
legend('Location', "southeast")


% example completed -------------------------------------------------------

completed = true;

% handling for ARCH competition
if nargout < 2
    clear res;
end
if nargout < 3
    clear tTotal;
end

end

% ------------------------------ END OF CODE ------------------------------
