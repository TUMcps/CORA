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
        [-0.1860741138812522; -0.0057845927755716; -0.0109108399869499; 0.0533508938505850; -0.3855531524843243], ...
        [0.3009754230893926; 0.4919841705605879; 0.1108221535648846; 0.2469338725012517; -0.1145206456580537]);;

% final result
assert(isequal(IH,IH_saved,1e-8));


% Reachability Analysis (zonotope bundles) --------------------------------

% redefine initial set
Z0{1} = params.R0;
Z0{2} = params.R0 + [-0.1; 0; 0.1; 0; 0];
params.R0 = zonoBundle(Z0);

R = reach(fiveDimSys, params, options);

IH = interval(R.timeInterval.set{end});

% saved result
IH_saved = interval( ...
        [-0.1860740281075652; -0.0057855253447022; -0.0109108315466627; 0.0533509231713888; -0.3855531524843243], ...
        [0.3005907856663744; 0.4913962991463204; 0.1108221535777970; 0.2469338725141641; -0.1145206456580537]);

% final result
assert(isequal(IH,IH_saved,1e-8));

res = true;

% ------------------------------ END OF CODE ------------------------------
