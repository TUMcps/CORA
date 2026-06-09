function res = testMosek_linearSysDT_observe_gain_06_PRadE()
% testMosek_linearSysDT_observe_gain_06_PRadE - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_06_PRadE
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Reference:
%    [1] Ye Wang, Meng Zhou, Vicenc Puig, Gabriela Cembrano, and
%        Zhenhua Wang. Zonotopic fault detection observer with H −
%        performance. In Proc. of the 36th IEEE Chinese Control
%        Conference, pages 7230–7235, 2017.

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
options.alg = 'PRad-E';

% observe
estSet = observe(vehicle,params,options);

% enclose last estimated set by interval
IH = interval(estSet.timePoint.set{end});

% Define comparison interval hull
IH_saved = interval([-0.217948002484438;-0.444592133619027;0.600833205200170;14.617475074845400;-0.222624181448318;-7.132392351168697],[0.188303672773522;-0.141467260451936;0.818454688867921;15.023020089590203;0.142908618751027;8.786563346538394]);

%check if slightly bloated versions enclose each other
assert(isequal(IH,IH_saved,1e-6));


end

% ------------------------------ END OF CODE ------------------------------
