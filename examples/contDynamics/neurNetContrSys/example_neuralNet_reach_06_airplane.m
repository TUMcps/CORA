function [completed,res,tTotal] = example_neuralNet_reach_06_airplane
% example_neuralNet_reach_06_airplane - example of reachability analysis
%    for a neural network controlled airplane                                  
%
% Syntax:
%    completed = example_neuralNet_reach_06_airplane()
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 
%    res - verification result
%    tTotal - total time
% 
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP21 Category Report: 
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants." 
%       EPiC Series in Computing 80 (2021): 90-119.

% Authors:       Niklas Kochdumper, Tobias Ladner
% Written:       08-November-2021
% Last update:   23-May-2022 (TL, ARCH'22 revisions)
%                30-March-2023 (TL, verify violated runs, ARCH'23 revisions)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp("BENCHMARK: Airplane")
rng(84)

% Parameters --------------------------------------------------------------

size = [0, 1]; % we show violation through simulation
R0 = [ % initial set
  0, 0; % x
  0, 0; % y
  0, 0; % z
  size; % u
  size; % v
  size; % w
  size; % phi
  size; % theta
  size; % psi
  0, 0; % r
  0, 0; % p
  0, 0; % q
];
R0 = interval(R0(:, 1), R0(:, 2));

params.tFinal = 2;
params.R0 = polyZonotope(R0);

% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.alg = 'lin';
options.tensorOrder = 2;
options.taylorTerms = 4;
options.zonotopeOrder = 50;

% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.poly_method = "singh";

% System Dynamics ---------------------------------------------------------

% open-loop system
sys = nonlinearSys(@airplane,12,6);

% load neural network controller
% [12, 100, 100, 20, 6]
nn = neuralNetwork.readONNXNetwork('controller_airplane.onnx');

% construct neural network controlled system
sys = neurNetContrSys(sys,nn,0.1);

% Specification -----------------------------------------------------------

% A1 = [0 1 0; 0 -1 0]; b1 = [0.5; 0.5];
% A2 = [eye(3); -eye(3)]; b2 = ones(6,1);
% set = polytope([blkdiag(A1,A2),zeros(8,6)],[b1;b2]);
safeSet = [
  -inf, inf; % x
  -0.5, 0.5; % y
  -inf, inf; % z
  -inf, inf; % u
  -inf, inf; % v
  -inf, inf; % w
  -1, 1; % phi
  -1, 1; % theta
  -1, 1; % psi
  -inf, inf; % r
  -inf, inf; % p
  -inf, inf; % q
];

safeSet = interval(safeSet(:, 1), safeSet(:, 2));
spec = specification(safeSet,'safeSet');

% Verification ------------------------------------------------------------

t = tic;
[res, R, simRes] = verify(sys, spec, params, options, evParams, true);
tTotal = toc(t);
disp(['Result: ' res])

% Visualization -----------------------------------------------------------

disp("Plotting..")
figure; hold on; box on;
projDims = [2,7];

% plot specification
plot(specification(safeSet, 'safeSet'), projDims, 'DisplayName', 'Safe set');

% plot reachable set
useCORAcolors("CORA:contDynamics")
plot(R, projDims, 'DisplayName', 'Reachable set')

% plot initial set
plot(R0, projDims, 'k', 'FaceColor', [1 1 1], 'DisplayName', 'Initial set');

% plot simulation
plot(simRes,projDims, 'DisplayName', 'Simulations');

% labels and legend
xlabel('y'); ylabel('\phi');
legend();


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
