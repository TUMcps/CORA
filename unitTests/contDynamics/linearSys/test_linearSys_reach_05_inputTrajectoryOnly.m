function res = test_linearSys_reach_05_inputTrajectoryOnly
% test_linearSys_reach_05_inputTrajectoryOnly - unit test for linear 
%    reachability analysis with an input trajectory uTransVec; this test 
%    should check whether correct the input trajectory is correctly
%    considered
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

% model parameters
params.tFinal = 0.12;
params.R0 = zonotope(ones(5,1));
params.u = 100*[1; 0; 0; 0.5; -0.5]; % start of input trajectory
params.u(:,2:3) = 0;
params.U = zonotope(zeros(5,2));

% reachability settings
options.timeStep = 0.04;
options.taylorTerms = 4;
options.zonotopeOrder = 200;

% system dynamics
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
fiveDimSys = linearSys('fiveDimSys',A,B);

% reachability analysis
R = reach(fiveDimSys, params, options);

% interval hull of final set
IH = interval(R.timeInterval.set{end});

% saved result
% IH_true = interval( ...
% [3.7035909068940835; 2.0587487720830868; 0.9218120355943664; 2.0792489149905693; -0.9222005200165104], ...
% [4.2545447422606566; 2.6163966074240355; 0.9494850143887323; 2.3899927521345345; -0.8505317877218803]);
IH_true = interval( ...
[3.7035968788199916; 2.0587560069906594; 0.9218129363717922; 2.0792497885965120; -0.9222004605735724], ...
[4.2545375074129845; 2.6163894787794950; 0.9494842613996752; 2.3899918745249180; -0.8505318311952287]);

% compare
assert(isequal(IH,IH_true,1e-8));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
