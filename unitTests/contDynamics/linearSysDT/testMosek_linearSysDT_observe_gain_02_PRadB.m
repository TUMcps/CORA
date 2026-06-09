function res = testMosek_linearSysDT_observe_gain_02_PRadB()
% testMosek_linearSysDT_observe_gain_02_PRadB - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_02_PRadB
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] V. T. H. Le, C. Stoica, T. Alamo, E. F. Camacho, and
%        D. Dumur. Zonotope-based set-membership estimation for
%        multi-output uncertain systems. In Proc. of the IEEE
%        International Symposium on Intelligent Control (ISIC),
%        pages 212–217, 2013.

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
options.alg = 'PRad-B';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.188666149875793;-0.247987569988081;0.546943766268291;14.694563460704801;-0.152604570119711;-7.648750940654575],[0.226554932974356;-0.230355071419466;0.750696505517899;14.954255177996770;0.110396172886972;11.514430420772320]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));

end

% ------------------------------ END OF CODE ------------------------------
