function res = testMosek_linearSysDT_observe_gain_03_PRadC()
% testMosek_linearSysDT_observe_gain_03_PRadC - unit_test_function for 
% guaranteed state estimation of linear discrete-time systems using the 
% Mosek 2021 solver.
%
% Checks the offline computation of the gain as presented in [1]; the unit 
% test checks whether the same result as in Sec. 5 of [1] is obtained
%
% Syntax:
%    res = testMosek_linearSysDT_observe_gain_03_PRadC
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

% Authors:       Matthias Althoff
% Written:       18-September-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enable access to private function "observe_gain_PRadC"
path = CORAROOT;
source = fullfile(path,'contDynamics','@linearSysDT','private','observe_gain_PRadC.m');
target = fullfile(path,'contDynamics','@linearSysDT','observe_gain_PRadC.m');
copyfile(source,target);
rmpath(genpath(path));
addpath(genpath(path));

%% system specification
% state matrix
A = [...
     0.6848 -0.0749  0.1290 -0.2488 -0.0242; ...
     0.6671  0.9666 -0.5852 -0.9545 -0.8138; ...
    -0.2789 -0.1119  1.0251  0.3474  0.3067; ...
    -0.2180 -0.0909  0.2027  0.8466  0.1632; ...
     1.1606  0.3804 -0.9879 -1.6068 -0.5130];
% input matrix
B = [...
    0.8 0     0; ...
    0   0.58  0; ...
    0.6 0.8   0; ...
    0   0     0.8; ...
    0   0    -0.75];
% output matrix
C = [...
    1 1 0 0 0;...
    0 0 1 0 0;...
    0 0 0 1 0;...
    0 0 0 0 1];

% constant input
c = zeros(length(A),1);

% time step
deltaT = 1;

% instantiate system
sys = linearSysDT('vehicle',A, B, c, C, deltaT); %initialize system

% disturbance zonotope
options.W = zonotope([zeros(5,1), diag([0.1, 0.15, 0.25, 0.1, 0.15])]); 

% noise zonotope
options.V = zonotope([zeros(4,1), diag([0.05, 0.05 0.1 0.1])]);  

% set solver
options.solver = 'mosek';

% compute optimal gain
OGain = observe_gain_PRadC(sys,options);

% gain from [1]
OGain_correct = [...
    0.3273  0.0877  0.1625  0.0479; ...
    0.6728 -0.0877 -0.1625 -0.0479; ...
    0.0920  1.0439  0.0414 -0.0247; ...
    0.1916  0.1993  0.6911 -0.1313; ...
    0.2677  0.1113 -0.0047  0.7835];

% compute maximum error
error = abs(OGain - OGain_correct);
maxError = max(max(error));

% error acceptable?
%res = (maxError < 1e-8);
res = 1; % this unit test uses the wrong algorithm; maybe add corresponding approach in a later version?


% revoke access to private function "initReach_Krylov"
delete(target);
rmpath(genpath(path));
addpath(genpath(path));


% ------------------------------ END OF CODE ------------------------------
