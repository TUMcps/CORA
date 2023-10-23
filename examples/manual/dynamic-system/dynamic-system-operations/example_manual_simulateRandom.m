function example_manual_simulateRandom()
% example_manual_simulateRandom - example from the manual demonstrating the 
% simulateRandom operation as defined in the manual
%
% Syntax:
%   example_manual_simulateRandom()
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

% system dynamics
sys = linearSys([-0.7 -2;2 -0.7],[1;1],[-2;-1]);

% parameter
params.tFinal = 5;
params.R0 = zonotope(interval([2;2],[2.5;2.5]));
params.U = zonotope(interval(-0.1,0.1));

% simulation settings
options.points = 7;
options.fracVert = 0.5;
options.fracInpVert = 1;
options.nrConstInp = 10;

% random simulation
simRes = simulateRandom(sys,params,options);

% reachability settings
options.timeStep = 0.05;
options.zonotopeOrder = 10;
options.taylorTerms = 5;

% reach
R = reach(sys,params,options);

% plot --------------------------------------------------------------------

figure; hold on;
useCORAcolors("CORA:contDynamics")
plot(R)
plot(R(1).R0)
plot(simRes)

% title('$\mathcal{S}$ and center','Interpreter','latex');
xlabel('$x_{(1)}$','Interpreter','latex')
ylabel('$x_{(2)}$','Interpreter','latex')

% ------------------------------ END OF CODE ------------------------------
