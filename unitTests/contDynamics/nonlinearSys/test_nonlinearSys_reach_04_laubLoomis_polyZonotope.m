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
 
% assume true
res = true;


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
    [1.0149071383850252; 0.8008249447366792; 0.9260067912066403; 1.6096886042592162; 0.4931654608095589; -0.0696709952280837; 0.0041738752897748], ...
    [1.6994210307009445; 1.5100585708305747; 1.6826517863619266; 2.4182337203368247; 1.1023506577076780; 0.2924850920995423; 0.7099585755176793]);

% final result
assert(isequal(IH,IH_saved,1e-8));

% ------------------------------ END OF CODE ------------------------------
