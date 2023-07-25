function completed = example_neuralNet_reach_06_airplane
% example_neuralNet_reach_06_airplane - example of reachability analysis
%                                       for a neural network controlled 
%                                       airplane                                  
%
% Syntax:  
%    completed = example_neuralNet_reach_06_airplane()
%
% Inputs:
%    no
%
% Outputs:
%    completed - true/false 
% 
% Reference:
%   [1] Johnson, Taylor T., et al. "ARCH-COMP21 Category Report: 
%       Artificial Intelligence and Neural Network Control Systems (AINNCS)
%       for Continuous and Hybrid Systems Plants." 
%       EPiC Series in Computing 80 (2021): 90-119.

% Author:       Niklas Kochdumper, Tobias Ladner
% Written:      08-November-2021
% Last update:  23-May-2022 (TL: ARCH'22 Revisions)
% Last revision:---

%------------- BEGIN CODE --------------

disp("BENCHMARK: Airplane")

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
options.taylorTerms = 4;
options.zonotopeOrder = 50;

options.alg = 'lin';
options.tensorOrder = 2;


% Parameters for NN evaluation --------------------------------------------

evParams = struct();
evParams.bound_approx = true;
evParams.polynomial_approx = "lin";


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
% set = mptPolytope([blkdiag(A1,A2),zeros(8,6)],[b1;b2]);
safeSet = [
  -inf, inf; % x
  -0.5, 0.5; % y
  -inf, inf; %
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


% Simulation --------------------------------------------------------------

tic
simRes = simulateRandom(sys, params);
tSim = toc;
disp(['Time to compute random simulations: ', num2str(tSim)]);


% Check Violation --------------------------------------------------------

tic
isVio = false;
for i = 1:length(simRes)
    x = simRes(i).x{1};
    for j=1:length(safeSet)
        isVio = isVio || ~all( ...
            (infimum(safeSet(j)) <= x(:, j)) & ...
            (x(:, j) <= supremum(safeSet(j))));
    end
end
tVio = toc;
disp(['Time to check violation in simulations: ', num2str(tVio)]);

if isVio
    disp("Result: VIOLATED")
    R = [];
    tComp = 0;
    tVeri = 0;
else
    % Reachability Analysis -----------------------------------------------

    tic
    R = reach(sys, params, options, evParams);
    tComp = toc;
    disp(['Time to compute reachable set: ', num2str(tComp)]);

    % Verification --------------------------------------------------------

    tic
    isVeri = true;
    for i = 1:length(R)
        R_i = R(i);
        for j = 1:length(R_i.timeInterval)
            isVeri = isVeri & safeSet.contains(R_i.timeInterval.set{j});
        end
    end
    tVeri = toc;
    disp(['Time to check verification: ', num2str(tVeri)]);

    if isVeri
        disp('Result: VERIFIED');
    else
        disp('Result: UNKNOWN');
    end
end

disp(['Total Time: ', num2str(tSim+tVio+tComp+tVeri)]);

% Visualization -----------------------------------------------------------

disp("Plotting..")
figure; hold on; box on;
projDims = [2,7];

% plot specification
plot(specification(safeSet, 'safeSet'), projDims, 'DisplayName', 'Safe set');

useCORAcolors("CORA:contDynamics")
if ~isVio
    % plot reachable set
    plot(R,projDims,'DisplayName','Reachable set');
    % plot initial set
    plot(R(1).R0,projDims,'DisplayName','Initial set');
else
    updateColorIndex()
    % plot initial set
    plot(R0,projDims,'DisplayName','Initial set', 'EdgeColor',[0 0 0],'FaceColor',CORAcolor("CORA:initialSet"));
end

% plot simulation
plot(simRes,projDims, 'DisplayName', 'Simulations');

% labels and legend
xlabel('y'); ylabel('\phi');
legend();


% example completed
completed = true;

%------------- END OF CODE --------------