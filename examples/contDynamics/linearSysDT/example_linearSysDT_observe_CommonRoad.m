function completed = example_linearSysDT_observe_CommonRoad
% example_linearSysDT_observe_CommonRoad - example for guaranteed state
%     estimation of linear discrete-time systems; the dynamics of traffic
%     participants from a CommonRoad file is estimated.
%
% Syntax:
%    completed = example_linearSysDT_observe_CommonRoad
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false

% Authors:       Matthias Althoff, Niklas Kochdumper, Carlos Valero
% Written:       15-July-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

%% System Dynamics

% continuous-time system
A = [0 0 1 0;0 0 0 1;0 0 0 0;0 0 0 0];
B = [0 0;0 0;1 0;0 1];
C = [1 0 0 0;0 1 0 0];

sys = linearSys(A,B,[],C);

% convert to discrete-time system
Ts = 0.01;                           % sampling time
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
if Ts ~= 0.1
    for i = 1:length(traj)
       tVec = 0:Ts:0.1*(size(traj{i}.x,2)-1);
       x = traj{i}.x;
       temp1 = interp1(traj{i}.t,x(1,:),tVec,'spline');
       temp2 = interp1(traj{i}.t,x(2,:),tVec,'spline');
       x = [temp1;temp2];
       traj{i}.t = tVec;
       traj{i}.x = x;
    end
end


%% Parameters

% measurement uncertainty
params.sigma = [0.1; 0.1];

% measurement 
params.y = traj{5}.x;

% initial set
dv = 50/3.6;
x0 = [traj{5}.x(:,1);0;0];
R0 = interval([-params.sigma;-dv;-dv],[params.sigma;dv;dv]) + x0;
params.R0 = zonotope(R0);

% set of model uncertainties
params.V = zonotope(zeros(2,1),0*eye(2));
params.W = zonotope(zeros(4,1),[zeros(2,2);11.5 0;0 11.5]);

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
plotOverTime(R,3);
xlabel('time');
ylabel('v_x');

figure; hold on;
plotOverTime(R,4);
xlabel('time');
ylabel('v_y');

figure; hold on;
plot(R,[1,2]);
plot(traj{5}.x(1,:),traj{5}.x(2,:),'Color',colorblind('y'));
xlabel('x');
ylabel('y');
axis equal;


% example completed
completed = true;


% ------------------------------ END OF CODE ------------------------------
