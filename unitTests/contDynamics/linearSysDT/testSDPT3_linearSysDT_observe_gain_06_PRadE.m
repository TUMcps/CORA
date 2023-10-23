function res = testSDPT3_linearSysDT_observe_gain_06_PRadE()
% testSDPT3_linearSysDT_observe_gain_06_PRadE - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testSDPT3_linearSysDT_observe_gain_06_PRadE
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "observe_gain_PRadE"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_PRadE.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_PRadE.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));
clear linearSysDT % required for some MATLAB versions so that copied file is recognized  

% load side slip model
load vehicleModel_dim2 vehicle params options

% exact result
PRadE = [ ...
0.0557145900132214; ...
0.0781710543629908];

% set solver
options.solver = 'sdpt3';

% compute optimal gain
options = params2options(params,options);
OGain = observe_gain_PRadE(vehicle,options);

% compute maximum error
error = abs(OGain - PRadE);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);


% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% ------------------------------ END OF CODE ------------------------------
