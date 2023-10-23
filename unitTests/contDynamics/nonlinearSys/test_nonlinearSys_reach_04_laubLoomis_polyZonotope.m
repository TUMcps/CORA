function res = test_nonlinearSys_reach_04_laubLoomis_polyZonotope
% test_nonlinearSys_reach_04_laubLoomis_polyZonotope - unit_test_function of 
%    nonlinear reachability analysis: Checks the solution of a 7D nonlinear
%    example using a non-convex set representation;
%
% Syntax:
%    res = test_nonlinearSys_reach_04_laubLoomis_polyZonotope
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       26-January-2016
% Last update:   23-April-2020 (restructure params/options)
%                22-June-2020 (NK, adapted to polyZonotope set rep.)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

x0 = [1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45];
R0 = zonotope([x0,0.3*eye(7)]);
params.R0 = polyZonotope(R0);                          % initial set
params.tFinal = 0.2;                                   % final time


% Reachability Settings ---------------------------------------------------

% settings
options.timeStep = 0.01;
options.taylorTerms = 20;
options.zonotopeOrder = 100;
options.intermediateOrder = 50;
options.errorOrder = 15;

% reachability algorithm
options.alg = 'poly';
options.tensorOrder = 3;

% settings for polynomial zonotopes
polyZono.maxDepGenOrder = 30;
polyZono.maxPolyZonoRatio = 0.001;
polyZono.restructureTechnique = 'reduceFullGirard';

options.polyZono = polyZono;


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@laubLoomis);


% Reachability Analysis ---------------------------------------------------

R = reach(sys, params, options);


% Numerical Evaluation ----------------------------------------------------

% enclose result by interval
IH = interval(R.timeInterval.set{end});

% saved result
IH_saved = interval( ...
    [1.014907148632731;0.800824936891000;0.926007014302413;1.609688394869659;0.493165360530853;-0.069670995159681;0.004173997644759], ...
    [1.699421020453211;1.510058578698442;1.682651563264945;2.418233929718181;1.102350758408782;0.292485092031139;0.709958453161779]);

% final result
res = isequal(IH,IH_saved,1e-8);

% ------------------------------ END OF CODE ------------------------------
