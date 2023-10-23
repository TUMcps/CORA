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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "observe_gain_HinfG"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_HinfG.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_HinfG.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));

% Load side slip model
load slipEstimationModel_6D params options vehicle

% Load exact result
load gain_sideSlip6D_Hinf-G HinfG

% compute optimal gain
options = params2options(params,options);
OGain = observe_gain_HinfG(vehicle,options);

% compute maximum error
error = abs(OGain - HinfG);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);


% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% ------------------------------ END OF CODE ------------------------------
