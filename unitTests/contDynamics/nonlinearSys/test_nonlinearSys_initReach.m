function res = test_nonlinearSys_initReach
% test_nonlinearSys_initReach - unit_test_function for computing a single
%    time interval reachable set for nonlinear dynamics:
%    Checks initReach of the nonlinearSys class for the 6 tank example;
%    It is checked whether partial reachable sets and the set
%    of linearization errors are correctly obtained
%
% Syntax:
%    res = test_nonlinearSys_initReach
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Matthias Althoff
% Written:       31-July-2017
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% model parameters
dim_x = 6;
params.R0 = zonotope([2; 4; 4; 2; 10; 4],0.2*eye(dim_x));
params.U = zonotope(0,0.005);
params.tFinal = 4;

% reachability settings
options.timeStep=4;
options.taylorTerms=4;
options.zonotopeOrder=50;
options.alg = 'lin';
options.tensorOrder = 2;

% system dynamics
tank = nonlinearSys(@tank6Eq);

% options check
[params,options] = validateOptions(tank,params,options,'FunctionName','reach');

% compute derivations (explicitly, since reach-function is not called)
derivatives(tank,options);

% obtain factors for reachability analysis
for i=1:(options.taylorTerms+1)
    %compute initial state factor
    options.factor(i) = options.timeStep^(i)/factorial(i);    
end

% comupute only first step
Rfirst = initReach(tank,params.R0,params,options);

% obtain interval hull of reachable set of first point in time
IH_tp = interval(Rfirst.tp{1}.set);
% obtain interval hull of reachable set of first time interval
IH_ti = interval(Rfirst.ti{1});
% obtain linearization errors
linErrors = Rfirst.tp{1}.error;


% ground truth
IH_tp_true = interval( ...
    [1.8057949711597598; 3.6433030183959114; 3.7940260617482671; 1.9519553317477598; 9.3409949650858550; 4.0928655724716370], ...
    [2.2288356782079028; 4.0572873081850807; 4.1960714210115002; 2.3451418924166987; 9.7630596270322201; 4.4862797486713282]);
IH_ti_true = interval( ...
    [1.7699801606999799; 3.6281401930838144; 3.7805292441504390; 1.7850641948695933; 9.3278848047801457; 3.7900869967590674], ...
    [2.2489804444442649; 4.2207006906857227; 4.2157311855735484; 2.3652362932387256; 10.2200546341070346; 4.5042192761237603]);
linErrors_true = 1e-3*[0.206863579523074; 0.314066666873806; 0.161658311464827; 0.353255431809860; 0.358487021465299; 0.209190642349808];

% compare results
assert(isequal(IH_tp,IH_tp_true,1e-8))
assert(isequal(IH_ti,IH_ti_true,1e-8))
assert(compareMatrices(linErrors,linErrors_true,1e-12));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
