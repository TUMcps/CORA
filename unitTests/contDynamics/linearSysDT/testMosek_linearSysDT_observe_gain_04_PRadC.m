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

% Author:       Matthias Althoff
% Written:      01-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% enable access to private function "observe_gain_PRadB"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_PRadC.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_PRadC.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));

% Load side slip model
load slipEstimationModel_6D params options vehicle

% Load exact results
load gain_sideSlip6D_PRad-C PRadC


%% PRadC
% compute optimal gain
options = params2options(params,options);
OGain = observe_gain_PRadC(vehicle,options);

% compute maximum error
error = abs(OGain - PRadC);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);

% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

%------------- END OF CODE --------------
