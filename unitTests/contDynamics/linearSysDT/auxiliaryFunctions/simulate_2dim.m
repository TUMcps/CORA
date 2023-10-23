function simulate_2dim
% simulate_2dim - simulates a two-dimensional system to disturbances and 
% sensor noise
%
% Syntax:
%    simulate_2dim
%
% Inputs:
%    -
%
% Outputs:
%    -
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
% Written:       04-June-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
savepath = [CORAROOT filesep 'unitTests' filesep 'contDynamics' filesep 'linearSysDT' filesep 'models'];

%% Settings

options.zonotopeOrder = 20; %zonotope order
options.timeStep = 1;


%% System Dynamics

% system matrices (we select delta = 1)
A = [0 -0.5;...
    1 1];
B = 1;
c = zeros(2,1);
C = [-2 1];

sys = linearSysDT('sys',A, B, c, C, options.timeStep); %initialize system

%% Parameters

params.tFinal = 20; %final time
params.R0 = zonotope(zeros(2,1),3*eye(2)); %initial set
params.V = 0.2*zonotope([0,1]); % sensor noise set
params.W = 0.02*[-6; 1]*zonotope([0,1]); % disturbance set
params.u = zeros(2,1); %input transition


%% Simulation Settings

options.points = 1;
options.p_conf = 0.999; % probability that sample of normal distribution within specified set
options.type = 'gaussian';

% simulate result assuming Gaussian distributions
simRes = simulateRandom(sys, params, options);

%% obtain output values
for i=1:length(simRes.t{1})
    % create measurement noise
    v = randPoint(params.V,'gaussian',options.p_conf);
    % obtain output value
    params.y(:,i) = C*simRes.x{1}(i,:)' + v;
end

save([savepath '/' 'twoDimSys_new'], 'sys', 'params', 'options', 'simRes');


% ------------------------------ END OF CODE ------------------------------
