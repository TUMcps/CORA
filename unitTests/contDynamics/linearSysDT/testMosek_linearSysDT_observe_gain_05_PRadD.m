function res = testMosek_linearSysDT_observe_gain_05_PRadD()
% testMosek_linearSysDT_observe_gain_05_PRadD - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_05_PRadD
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] Ye Wang, Zhenhua Wang, Vicenc Puig, and Gabriela
%        Cembrano. Zonotopic set-membership state estimation for
%        discrete-time descriptor LPV systems. IEEE Transactions
%        on Automatic Control, 64(5):2092-2099, 2019.

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
options.alg = 'PRad-D';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.191856670961832;-0.245816042378656;0.528557553823015;14.717327962524285;-0.120735675623385;-13.659914995675745],[0.186191974746718;-0.232458222003942;0.728748494647000;14.919469699762045;0.082756833359093;14.131792763257149]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));

end

% ------------------------------ END OF CODE ------------------------------
