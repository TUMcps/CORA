function res = testMosek_linearSysDT_observe_gain_04_PRadC()
% testMosek_linearSysDT_observe_gain_04_PRadC - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1],[2]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_04_PRadC
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] Ye Wang, Teodoro Alamo, Vicenc Puig, and Gabriela
%        Cembrano. A distributed set-membership approach based on
%        zonotopes for interconnected systems. In Proc. of the IEEE
%        Conference on Decision and Control (CDC), pages 668–673, 2018.
%    [2] Ye Wang, Vicenç Puig, and Gabriela Cembrano. Set-
%        membership approach and kalman observer based on
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
options.alg = 'PRad-C';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.184146962945192;-0.247266897577215;0.526779736482237;14.717179472433656;-0.120534778597817;-9.177810294598721],[0.195604585556369;-0.231133396380272;0.727075472671465;14.918706347375380;0.082214161094085;10.955914749279398]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));

end

% ------------------------------ END OF CODE ------------------------------
