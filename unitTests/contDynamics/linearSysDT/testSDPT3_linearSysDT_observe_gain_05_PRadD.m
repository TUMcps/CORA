function res = testSDPT3_linearSysDT_observe_gain_05_PRadD()
% testSDPT3_linearSysDT_observe_gain_05_PRadD - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:  
%    res = testSDPT3_linearSysDT_observe_gain_05_PRadD
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

% Author:       Matthias Althoff
% Written:      01-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% enable access to private function "observe_gain_PRadD"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_PRadD.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_PRadD.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));
clear linearSysDT % required for some MATLAB versions so that copied file is recognized  

% load side slip model
load vehicleModel_dim2 vehicle params options

% exact result
PRadD = [ ...
0.0210427394634470; ...
0.0666743436391924];

% set solver
options.solver = 'sdpt3';

% compute optimal gain
options = params2options(params,options);
OGain = observe_gain_PRadD(vehicle,options);

% compute maximum error
error = abs(OGain - PRadD);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);


% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

%------------- END OF CODE --------------
