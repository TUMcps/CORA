function completed = example_linearSysDT_observe_02_2dim
% example_linearSysDT_observe_02_2dim - example for guaranteed state
%    estimation of linear discrete-time systems from a unit test; shows the
%    solution of the linearSysDT class for a two-dimensional example
%    from Sec. 7.1 of [1].
%
% Syntax:
%    res = example_linearSysDT_observe_02_2dim
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043, 2005.

% Authors:       Matthias Althoff
% Written:       30-April-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 20; %final time
params.R0 = zonotope(zeros(2,1),3*eye(2)); %initial set
params.V = 0.2*zonotope([0,1]); % sensor noise set
params.W = 0.02*[-6; 1]*zonotope([0,1]); % disturbance set
params.u = zeros(2,20); % input vector
params.y = [0.79, 5.00, 4.35, 1.86, -0.11, -1.13, -1.17, -0.76, ...
    -0.12, 0.72, 0.29, 0.19, 0.09, -0.21, 0.05, -0.00, -0.16, 0.01, ...
    -0.08, 0.13]; %measurement vector


% Algorithmic Settings ----------------------------------------------------

options.zonotopeOrder = 20; % zonotope order
options.timeStep = 1; % setp size
options.alg = 'FRad-C'; % observer approach


% System Dynamics ---------------------------------------------------------

A = [0 -0.5; 1 1];
B = 1;
c = zeros(2,1);
C = [-2 1];
twoDimSys = linearSysDT('twoDimSys',A, B, c, C, options.timeStep); 


% Set-based observation ---------------------------------------------------

EstSet = observe(twoDimSys,params,options);


% Visualization -----------------------------------------------------------

for iDim = 1:2
    figure; hold on;
    % plot time elapse
    plotOverTime(EstSet,iDim);

    % label plot
    xlabel('t');
    ylabel(['x_{',num2str(iDim),'}']);
end

completed = true;

% ------------------------------ END OF CODE ------------------------------
