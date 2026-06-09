function res = testMosek_linearSysDT_observe_gain_08_HinfG()
% testMosek_linearSysDT_observe_gain_08_HinfG - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_08_HinfG
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] W. Tang, Z. Wang, Y. Wang, T. Raissi, and Y. Shen.
%        Interval estimation methods for discrete-time linear time-
%        invariant systems. IEEE Transactions on Automatic Control,
%        64(11):4717-4724, 2019.

% Authors:       Matthias Althoff
% Written:       01-March-2021
% Last update:   15-April-2026 (NH, changed to observe(sys,params,options))
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% Load model
[vehicle,params,options] = load_model_linearSysDT("slipEstimationModel_6D");

% Select observer
options.alg = 'Hinf-G';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.234655075321623;-0.481579843251740;0.554658876207267;14.588346676484386;-0.259761496741450;-3.659319275873318],[0.154502314705352;-0.170185044057171;0.788652846264336;15.059523221350046;0.134626789899628;3.936102530810363]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));

end

% ------------------------------ END OF CODE ------------------------------
