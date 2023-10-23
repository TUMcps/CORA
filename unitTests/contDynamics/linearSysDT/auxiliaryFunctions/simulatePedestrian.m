function simulatePedestrian
% simulatePedestrian - simulates a pedestrian subject to disturbances and 
% sensor noise
%
% Syntax:
%    simulatePedestrian(options)
%
% Inputs:
%    ---
%
% Outputs:
%    sys - system struct
%    options - options struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Authors:       Matthias Althoff
% Written:       03-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
savepath = [CORAROOT filesep 'unitTests' filesep 'contDynamics' filesep 'linearSysDT' filesep 'models'];

%% Settings

options.timeStep = 0.01;
options.zonotopeOrder = 100;
options.reductionTechnique = 'pca';

%% Input trajectory 
% We use jerk templates for changing acceleration in x- and y-direction

% input patterns
segmentDuration = 0.1; %[s]
segmentSteps = ceil(segmentDuration/options.timeStep);
jerk = 10/3; %[m/s^3]
noJerk = zeros(2,segmentSteps); % more jerk in x-direction
incrJerk_x = jerk*[ones(1,segmentSteps); zeros(1,segmentSteps)]; % positive jerk in x-direction
decrJerk_x = -jerk*[ones(1,segmentSteps); zeros(1,segmentSteps)]; % negative jerk in x-direction
incrJerk_y = jerk*[zeros(1,segmentSteps); ones(1,segmentSteps)]; % positive jerk in y-direction
decrJerk_y = -jerk*[zeros(1,segmentSteps); ones(1,segmentSteps)]; % negative jerk in y-direction

% jerk
jerk = [incrJerk_x, incrJerk_x, incrJerk_x, incrJerk_y, decrJerk_x, ...
    incrJerk_y, decrJerk_x, decrJerk_x, decrJerk_y, decrJerk_y, noJerk, ...
    noJerk, noJerk, noJerk, noJerk, noJerk, noJerk, noJerk, decrJerk_x, ...
    decrJerk_y, decrJerk_x, decrJerk_y, decrJerk_x, decrJerk_y];

% accerleration
acceleration = [0;0]; 
for i = 1:length(jerk(1,:))
    acceleration(:,i+1) = acceleration(:,i) + options.timeStep*jerk(:,i);
end

%% System Dynamics 
% system matrices
A = [ ...
0, 0, 1, 0; ...
0, 0, 0, 1; ...
0, 0, 0, 0; ...
0, 0, 0, 0];

B = [ ...
0, 0; ...
0, 0; ...
1, 0;...
0, 1];

A_d = expm(A*options.timeStep);

B_d = (eye(4)*options.timeStep + A*options.timeStep^2/2)*B;

C = [1 0 0 0; 0 1 0 0];

c = zeros(length(A), 1);


pedestrian = linearSysDT('pedestrian',A_d, B_d, c, C, options.timeStep); %initialize system


%% Parameters

% create unit ball for a_max<1 [m/s^2]
E_ball = ellipsoid(diag([0, 0, 1, 1]));
Z_ball = zonotope(E_ball, 10, 'outer:norm');
W = Z_ball; % disturbance set
W_d = (eye(4)*options.timeStep + A*options.timeStep^2/2)*W; % discrete-time disturbance set
V = zonotope([zeros(2,1), 0.1*eye(2)]); % sensor noise set

params.tFinal = length(jerk(1,:))*options.timeStep; %final time
params.R0 = zonotope(zeros(4,1),eye(4)); %initial set
params.V = V; % sensor noise set
params.W = W_d; % disturbance set
params.u = acceleration; %input transition


%% Simulation Settings

options.points = 1;
options.p_conf = 0.999; % probability that sample of normal distribution within specified set
options.type = 'gaussian';

% simulate result assuming Gaussian distributions
simRes = simulateRandom(pedestrian, params, options);

%% obtain output values
for i=1:length(simRes.t{1})
    % create measurement noise
    v = randPoint(params.V,'gaussian',options.p_conf);
    % obtain output value
    params.y(:,i) = C*simRes.x{1}(i,:)' + v;
end

save([savepath '/' 'pedestrianModel_new'], 'pedestrian', 'params', 'options', 'simRes');

% ------------------------------ END OF CODE ------------------------------
