function res = test_linearSys_reach_02_5dim()
% test_linearSys_reach_02_5dim - unit test function of linear reachability 
%    analysis with uncertain inputs
%
% Checks the solution of the linearSys class for a 5-dimensional example;
% It is checked whether the enclosing interval of the final reachable set 
% is close to an interval provided by a previous solution that has been saved
%
% Syntax:
%    res = test_linearSys_reach_02_5dim()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       09-August-2016
% Last update:   23-April-2020 (restructure params/options)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

dim_x = 5;
params.tFinal = 5;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = zonotope([1; 0; 0; 0.5; -0.5],0.5*diag([0.2, 0.5, 0.2, 0.5, 0.5]));


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.04; % time step size for reachable set computation
options.taylorTerms = 4; % number of taylor terms in exponential matrix
options.zonotopeOrder = 10; % zonotope order
options.linAlg = 'wrapping-free'; % algorithm


% System Dynamics ---------------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;
fiveDimSys=linearSys('fiveDimSys',A,B);


% Reachability Analysis (zonotope) ----------------------------------------

R = reach(fiveDimSys, params, options);
IH = interval(R.timeInterval.set{end});


% saved result
IH_saved = interval( ...
           [-0.186078149456309; -0.004994826957236; -0.010811262167507; 0.053366186432215; -0.385353029993981], ...
           [0.300725848257715; 0.491957694383420; 0.110810877781333; 0.246634559097104; -0.114528793200091]);

% final result
res_zono = isequal(IH,IH_saved,1e-8);


% Reachability Analysis (zonotope bundles) --------------------------------

% redefine initial set
Z0{1} = params.R0;
Z0{2} = params.R0 + [-0.1; 0; 0.1; 0; 0];
params.R0 = zonoBundle(Z0);

R = reach(fiveDimSys, params, options);

IH = interval(R.timeInterval.set{end});

% saved result
IH_saved = interval( ...
    [-0.1860781494563091; -0.0049948269572356; -0.0108112537143071; 0.0533662157659313; -0.3853530299939805], ...
    [0.3003413164189264; 0.4913712245778811; 0.1108108777813331; 0.2466345590971040; -0.1145287932000914]);

% final result
res_zonoBundles = isequal(IH,IH_saved,1e-8);

% result of different set representations
res = res_zono && res_zonoBundles;

% ------------------------------ END OF CODE ------------------------------
