function res = testMosek_linearSysDT_observe_gain_01_PRadA()
% testMosek_linearSysDT_observe_gain_01_PRadA - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in Sec. 4.1 of 
% [1]; the unit test checks whether the same result as in a previous
% implementation is obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_01_PRadA()
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] V. T. H. Le, C. Stoica, T. Alamo, E. F. Camacho, and
%        D. Dumur. Zonotopic guaranteed state estimation for
%        uncertain systems. Automatica, 49(11):3418–3424, 2013.

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   15-April-2026 (NH, changed to observe(sys,params,options))
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% Load model
[vehicle,params,options] = load_model_linearSysDT("slipEstimationModel_6D");

% Select observer
options.alg = 'PRad-A';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.778056174558984;-1.326674879561634;0.200143289514439;12.473975086385328;-0.384359504384016;-16.441513231174888],[1.183980022309567;1.083577489492029;0.916131742263717;16.129019752785275;0.543922444305893;15.918231613685712]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));

end

% ------------------------------ END OF CODE ------------------------------
