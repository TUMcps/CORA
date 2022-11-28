function res = testSDPT3_linearSysDT_observe_gain_01_PRadA()
% testSDPT3_linearSysDT_observe_gain_01_PRadA - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the offline computation of the gain as presented in Sec. 4.1 of 
% [1]; the unit test checks whether the same result as in a previous
% implementation is obtained
%
% Syntax:  
%    res = testSDPT3_linearSysDT_observe_gain_01_PRadA()
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

% Author:       Matthias Althoff
% Written:      18-September-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% enable access to private function "observe_gain_PRadA"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_PRadA.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_PRadA.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));
clear linearSysDT % required for some MATLAB versions so that copied file is recognized  

% load side slip model
load vehicleModel_dim2 vehicle params options

% exact result
PRadA = [ ...
0.0566853336654368; ...
0.0664455001897976];

% set solver
options.solver = 'sdpt3';

% compute optimal gain
options = params2options(params,options);
OGain = observe_gain_PRadA(vehicle,options);

% compute maximum error
error = abs(OGain - PRadA);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);


% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

%------------- END OF CODE --------------
