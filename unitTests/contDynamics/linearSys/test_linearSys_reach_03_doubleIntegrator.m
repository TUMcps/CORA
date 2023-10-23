function res = test_linearSys_reach_03_doubleIntegrator()
% test_linearSys_reach_03_doubleIntegrator - unit test for linear reachability 
% analysis with uncertain inputs; this test should check whether correct
% results are returned when the system matrix only consists of zeros
%
% Syntax:
%    res = test_linearSys_reach_03_doubleIntegrator()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       12-November-2017
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% System Dynamics ---------------------------------------------------------

A = zeros(2);
B = 1;
doubleIntegrator = linearSys('twoDimSys',A,B); %initialize system
dim_x = length(A);


% Parameters --------------------------------------------------------------

params.tFinal = 1;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = zonotope([1; 1],diag([0.1, 0.1]));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.04; % time step size for reachable set computation
options.taylorTerms = 4; % number of taylor terms for exponential matrix
options.zonotopeOrder = 10; % zonotope order
options.linAlg = 'wrapping-free';


% Reachability Analysis (zonotope) ----------------------------------------

Rcont = reach(doubleIntegrator, params, options);

IH = interval(Rcont.timeInterval.set{end});

% saved result
IH_true = ones(dim_x,1)*interval(1.76, 2.2);

% final result
res = isequal(IH,IH_true,1e-8);

% ------------------------------ END OF CODE ------------------------------
