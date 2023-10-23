function example_manual_observe()
% example_manual_observe - example from the manual demonstrating the 
% observe operation as defined in the manual
%
% Syntax:
%   example_manual_observe()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also:

% Authors:       Tobias Ladner
% Written:       27-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% Parameters
params.tFinal = 20; %final time
params.R0 = zonotope(zeros(2,1),3*eye(2)); %initial set
params.V = 0.2*zonotope([0,1]); % sensor noise set
params.W = 0.02*[-6; 1]*zonotope([0,1]); % disturbance set
params.u = zeros(2,1); % input vector
params.y = [0.79, 5.00, 4.35, 1.86, -0.11, -1.13, -1.17, -0.76, ...
    -0.12, 0.72, 0.29, 0.19, 0.09, -0.21, 0.05, -0.00, -0.16, 0.01, ...
    -0.08, 0.13]; %measurement vector

%% Algorithmic Settings
options.zonotopeOrder = 20; % zonotope order
options.timeStep = 1; % setp size
options.alg = 'FRad-C'; % observer approach

%% System Dynamics
reactor = linearSysDT('reactor',[0 -0.5; 1 1], 1, zeros(2,1), [-2 1], options.timeStep);

% observe
EstSet = observe(reactor,params,options);

% plot --------------------------------------------------------------------

figure;
subplot(1, 2, 1); hold on;
useCORAcolors("CORA:contDynamics")
plotOverTime(EstSet,1)

xlabel('time','Interpreter','latex')
ylabel('$x_{(1)}$','Interpreter','latex')

subplot(1, 2, 2); hold on;
useCORAcolors("CORA:contDynamics")
plotOverTime(EstSet,2)

xlabel('time','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
