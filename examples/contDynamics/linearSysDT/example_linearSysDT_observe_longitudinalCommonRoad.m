function completed = example_linearSysDT_observe_longitudinalCommonRoad
% example_linearSysDT_observe_longitudinalCommonRoad - example for 
%     guaranteed state estimation of linear discrete-time systems.
%     The longitundinal dynamics of traffic participants from a CommonRoad
%     file is estimated.
%
% Syntax:  
%    completed = example_linearSysDT_observe_longitudinalCommonRoad
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false 

% Author:       Matthias Althoff, Niklas Kochdumper, Carlos Valero
% Written:      15-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


%% System Dynamics

% continuous-time system
A = [0 1; 0 0];
B = [0;1];
C = [1 0];

sys = linearSys(A,B,[],C);

% convert to discrete-time system
Ts = 0.1;                           % sampling time
sys = linearSysDT(sys,Ts);



%% Load CommonRoad Scenario 
[~,dynObs] = commonroad2cora('ARG_Carcarana-10_1_T-1', ...
                                          'verbose', false);
                                      
% extract trjectories for the vehicles
traj = {};
counter = 1;
x = []; t = [];
while counter <= length(dynObs)
    if counter ~= 1 && infimum(dynObs{counter}.time) == 0 
        traj{end+1}.x = x;
        traj{end}.t = t; 
        x = []; t = [];
    end
    x = [x,center(dynObs{counter}.set)];
    t = [t,supremum(dynObs{counter}.time)];
    counter = counter + 1;
end
traj{end+1}.x = x;
traj{end}.t = t;

% interpolation to increase observer frequency
for i = 1:length(traj)
   tVec = 0:Ts:0.1*(size(traj{i}.x,2)-1);
   x = traj{i}.x;
   temp1 = interp1(traj{i}.t,x(1,:),tVec,'spline');
   temp2 = interp1(traj{i}.t,x(2,:),tVec,'spline');
   x = [temp1;temp2];
   traj{i}.t = tVec;
   traj{i}.x = temp1;
end



%% Parameters

% measurement uncertainty
params.sigma = 0.1;

% measurement 
params.y = traj{7}.x;

% initial set
dv = 100/3.6;
x0 = [traj{7}.x(:,1);0];
R0 = interval([-params.sigma;-dv],[params.sigma;dv]) + x0;
params.R0 = zonotope(R0);

% set of model uncertainties
params.V = zonotope(zeros(1,1),0*eye(1));
params.W = zonotope([0,0; 0,11.5]);

% final time
params.tFinal = size(params.y,2)*Ts;


%% Settings
options.zonotopeOrder = 10;
options.timeStep = Ts;
options.alg = 'FRad-C';


%% Observation
R = observe(sys,params,options);


%% Visualization
figure; hold on;
plotOverTime(R,2);
xlabel('time');
ylabel('v_x');


% example completed
completed = true;


%------------- END OF CODE --------------
