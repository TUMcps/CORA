function res = test_nonlinearSys_reach_02_linearEqualsNonlinear
% test_nonlinearSys_reach_02_linearEqualsNonlinear - unit_test_function for 
%     nonlinear reachability analysis; it is checked whether the solution
%     of a linear system equals the solution obtained by the nonlinear
%     system class when the system is linear
%
% Syntax:
%    res = test_nonlinearSys_reach_02_linearEqualsNonlinear
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
params.tFinal = 0.04;
params.R0 = zonotope(ones(dim_x,1),0.1*eye(dim_x));
params.U = zonotope([1; 0; 0; 0.5; -0.5],0.5*diag([0.2, 0.5, 0.2, 0.5, 0.5]));


% Reachability Settings (linear) ------------------------------------------

optionsLin.timeStep = 0.04;
optionsLin.taylorTerms = 8;
optionsLin.zonotopeOrder = 10;
optionsLin.linAlg = 'wrapping-free';


% Reachability Settings (nonlinear) ---------------------------------------

optionsNonLin.timeStep = optionsLin.timeStep;
optionsNonLin.taylorTerms = optionsLin.taylorTerms;
optionsNonLin.zonotopeOrder = optionsLin.zonotopeOrder;


% System Dynamics ---------------------------------------------------------

% linear system
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2]; B = 1;
fiveDimSys = linearSys('fiveDimSys',A,B);

% nonlinear system
fiveDimSysNonlinear = nonlinearSys(@fiveDimSysEq);


% Reachability Analysis ---------------------------------------------------

% linear system
Rlin = reach(fiveDimSys, params, optionsLin);

% nonlinear system (conservative linearization)
optionsNonLin.alg = 'lin';
optionsNonLin.tensorOrder = 2;

Rnonlin1 = reach(fiveDimSysNonlinear, params, optionsNonLin);

% nonlinear system (conservative polynomialization)
optionsNonLin.alg = 'poly';
optionsNonLin.tensorOrder = 3;
optionsNonLin.intermediateOrder = 100*optionsNonLin.zonotopeOrder;
optionsNonLin.errorOrder = 1;

Rnonlin2 = reach(fiveDimSysNonlinear, params, optionsNonLin);


% Numerical Evaluation ----------------------------------------------------

% enclosing intervals of final reachable sets
IH = interval(Rlin.timePoint.set{end});
IH_nonlinear_T1 = interval(Rnonlin1.timePoint.set{end});
IH_nonlinear_T2 = interval(Rnonlin2.timePoint.set{end});

% final result
res = isequal(IH,IH_nonlinear_T1,1e-8) && isequal(IH,IH_nonlinear_T2,1e-8);

% ------------------------------ END OF CODE ------------------------------
