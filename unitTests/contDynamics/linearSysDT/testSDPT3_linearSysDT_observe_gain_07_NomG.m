function res = testSDPT3_linearSysDT_observe_gain_07_NomG()
% testSDPT3_linearSysDT_observe_gain_07_NomG - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% SDPT3 solver (version 4.0).
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in a previous implementation is 
% obtained
%
% Syntax:
%    res = testSDPT3_linearSysDT_observe_gain_07_NomG
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "observe_gain_NomG"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_NomG.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_NomG.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));
clear linearSysDT % required for some MATLAB versions so that copied file is recognized  

% load side slip model
load vehicleModel_dim2 vehicle params options

% exact result
NomG = [ ...
0.0165139190565368; ...
0.0646584371768353];

% set solver
options.solver = 'sdpt3';

% compute optimal gain
OGain = observe_gain_NomG(vehicle,options);

% compute maximum error
error = abs(OGain - NomG);
maxError = max(max(error));

% error acceptable?
res = (maxError < 1e-8);


% revoke access to private function
delete(target);
rmpath(genpath(path));
addpath(genpath(path));

% ------------------------------ END OF CODE ------------------------------
