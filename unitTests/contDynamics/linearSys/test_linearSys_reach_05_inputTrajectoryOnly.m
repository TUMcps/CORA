function res = test_linearSys_reach_05_inputTrajectoryOnly()
% test_linearSys_reach_05_inputTrajectoryOnly - unit test for linear 
% reachability analysis with an input trajectory uTransVec; this test 
% should check whether correct the input trajectory is correctly considered
%
% Syntax:
%    res = test_linearSys_reach_05_inputTrajectoryOnly
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       18-September-2018
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------
dim_x = 5;
params.tFinal = 0.12; %final time
params.R0 = zonotope(ones(dim_x,1),zeros(dim_x,1));
params.u = 100*[1; 0; 0; 0.5; -0.5]; % start of input trajectory
params.u(:,2:3) = 0;
params.U = zonotope(zeros(5,2));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.04; % time step size for reachable set computation
options.taylorTerms = 4; % number of taylor terms for exponential matrix
options.zonotopeOrder = 200; % zonotope order


% System Dynamics ---------------------------------------------------------
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
fiveDimSys = linearSys('fiveDimSys',A,B);


% Reachability Analysis (zonotope) ----------------------------------------
Rcont = reach(fiveDimSys, params, options);

% interval hull of final set
IH = interval(Rcont.timeInterval.set{end});

% saved result
IH_true = interval( ...
[3.7035909068940835; 2.0587487720830868; 0.9218120355943664; 2.0792489149905693; -0.9222005200165104], ...
[4.2545447422606566; 2.6163966074240355; 0.9494850143887323; 2.3899927521345345; -0.8505317877218803]);

% final result
res = isequal(IH,IH_true,1e-8);

% ------------------------------ END OF CODE ------------------------------
