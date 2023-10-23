function completed = example_linearSysDT_reach_07_10dim_backOver()
% example_linearSysDT_reach_07_10dim_backOver - example of discrete-time linear 
%    backward reachability analysis with uncertain inputs taken from 
%    Sec. VI.B in [1].
%
% Syntax:
%    example_linearSysDT_reach_07_10dim_backOver
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] L. Yang and N. Ozay, "Scalable Zonotopic Under-Approximation of 
%        Backward Reachable Sets for Uncertain Linear Systems," in IEEE 
%        Control Systems Letters, vol. 6, pp. 1555-1560, 2022.

% Authors:       Matthias Althoff
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Warning
disp('Warning: This example typically runs for many hours. Please comment the following lines to run it.');
completed = false;
return

% System Dynamics ---------------------------------------------------------

% system matrix of continuous dynamics
A_cont = ...
    [0 1 0 0 0 0 1     0  0 1; ... % 
    0 0 0 0 0 0 0     0  0 0; ...
    0 0 0 1 0 0 0    -1  0 0; ...
    0 0 0 0 0 0 0     0  0 0; ...
    0 0 0 0 0 1 0     0  1 0; ...
    0 0 0 0 0 0 0     0  0 0; 
    0 0 0 0 0 1 -0.01 1  0 0; ...
    0 0 0 0 0 0 -0.01 -1 0 0; ...
    0 0 0 0 0 0 -1e-4 0  0 2; ...
    0 0 0 0 0 0 0     -2 0 -1e-4];

% sampling time
dt = 0.5;

% system matrix
A = expm(A_cont*dt);

% input matrix of continuous dynamics
B_cont = ...
    [0 0 0;...
    1 0 0;...
    0 0 0;...
    0 1 0;...
    0 0 0;...
    0 0 1;...
    0 0 0;... 
    0 0 0;... 
    0 0 0;...
    0 0 0]; 

% input matrix approximated by 20th order Taylor approximation of the 
% integral of e^{At}; required becuase A is not invertible
A_int = zeros(length(A)); % initial value
for i=1:20
    A_int = A_int + A^(i-1)*dt^i/factorial(i);
end
B = A_int*B_cont;

% constant input
c = zeros(length(A),1);

% instantiate linear discrete time system
sys = linearSysDT('sys',A,B,c,dt);


% Parameter ---------------------------------------------------------------

params.tFinal = 25;
params.R0 = zonotope([[10; 0; 10; 0; 10; zeros(5,1)],0.5*eye(length(A))]);
params.U = zonotope([zeros(3,1), 0.5*eye(3)]);
params.W = zonotope([zeros(10,1),diag([0.12,0.2,0.12,0.2,0.12,0.2,0.1*ones(1,4)])]);


% Reachability Settings ---------------------------------------------------

%options.zonotopeOrder = 20;
options.zonotopeOrder = 15;
options.linAlg = 'backward_maxmin_coarse';
%options.linAlg = 'backward_maxmin_RKweighted';


% Reachability Analysis ---------------------------------------------------

tic
R = reach(sys, params, options);
tComp = toc;
disp(['computation time of reachable set: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

% final reachable set becomes new initial set
params.R0 = R.timePoint.set{end};
simOpt.points = 25;
simOpt.type = 'constrained';
simOpt.R = R;
simRes = simulateRandom(sys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1 2]};

for k = 1:length(dims)
    
    figure; hold on; box on
    useCORAcolors("CORA:contDynamics")
    projDims = dims{k};

    % plot reachable set
    plot(R,projDims);
    
    % plot initial output set
    plot(R.R0,projDims);
    
    % plot simulation results
    plot(simRes,projDims,'Marker','.');

    % label plot
    xlabel(['x_{',num2str(projDims(1)),'}']);
    ylabel(['x_{',num2str(projDims(2)),'}']);
end

% example completed
completed = true;

% ------------------------------ END OF CODE ------------------------------
