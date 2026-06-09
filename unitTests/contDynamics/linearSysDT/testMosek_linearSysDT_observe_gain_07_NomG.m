function res = testMosek_linearSysDT_observe_gain_07_NomG()
% testMosek_linearSysDT_observe_gain_07_NomG - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_07_NomG
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and Kalman observer based on
%        zonotopes for discrete-time descriptor systems. Automatica,
%        93:435-443, 2018.

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
options.alg = 'Nom-G';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.225898704526925;-0.443530742956009;0.627544324295057;14.656518527163982;-0.217104466264740;-9.547719910992430],[0.210964292493496;-0.127758196118713;0.838615415997879;15.014442017597094;0.187817674437501;10.464502863326619]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));

end

% ------------------------------ END OF CODE ------------------------------
