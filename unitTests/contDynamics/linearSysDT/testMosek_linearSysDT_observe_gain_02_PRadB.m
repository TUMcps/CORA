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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "observe_gain_PRadB"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_PRadB.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_PRadB.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));

% Load side slip model
load slipEstimationModel_6D params options vehicle

% Load exact result
load gain_sideSlip6D_PRad-B PRadB

% compute optimal gain
options = params2options(params,options);
OGain = observe_gain_PRadB(vehicle,options);

% compute maximum error
error = abs(OGain - PRadB);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);


% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% ------------------------------ END OF CODE ------------------------------
