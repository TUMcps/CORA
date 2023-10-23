function simulateVehicle
% simulateVehicle - simulates a vehicle subject to disturbances and sensor
% noise
%
% Syntax:
%    simulateVehicle
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
% Written:       08-September-2020
% Last update:   10-January-2021
%                04-June-2021
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set path
savepath = [CORAROOT filesep 'unitTests' filesep 'contDynamics' filesep 'linearSysDT' filesep 'models'];

% model dimension
dim = 6;

%% Settings

options.timeStep = 0.01;
options.zonotopeOrder = 100;
options.reductionTechnique = 'pca';

%% Input trajectory

% Steering magnitude
Nmag1 = 0.1; 

% input signal for double lane change
InpSig = [zeros(50,1);Nmag1*ones(50,1);-Nmag1*ones(50,1); zeros(50,1); ...
    -Nmag1*ones(50,1); Nmag1*ones(50,1); zeros(50,1)];

%% System Dynamics

% assign nr of outputs to different vehicle models 
if dim == 2
    nrOfOutputs = 1;
elseif dim == 4
    nrOfOutputs =3;
elseif dim == 6
    nrOfOutputs =4;
end

%set sensor noise
Vk = 0.1*ones(1,6);
% set disturbance
Wk = [0.005,0.055,0.005,0.04,0.02,0.2];

% vehicle paramters from CommonRoad vehicle paramters considering mu = 1
mu=1; 
Vx = 15;
lr = 1.508; 
lf = 0.883; 
m =1.225e3 ; 
g = 9.81;  
Iz = 1.538e3;
Cf = 45000; 
Cr = 57000;  
Bs = 0.5; 
Js = 0.05; 
etat = 0.010; 
Rs = 16;
Fzf = m*g*lr/(lf+lr); % front tire vertical force
Fzr = m*g*lf/(lf+lr); % rear tire vertical force
Csf = Cf/(mu*Fzf);Csr =Cr/(mu*Fzr);

% system dynamics
C1 = mu/(Vx*(lr+lf)); 
C2 = mu*m/(Iz*(lr+lf)); 
VehParam.a11 = -C1*g*(Csr*lf-Csf*lr);
VehParam.a12 = C1*(g*lf*lr/Vx)*(Csr-Csf) - 1; 
VehParam.a21 = C2*g*lr*lf*(Csr-Csf); 
VehParam.a22 = -C2*(g/Vx)*(lf*lf*lr*Csf + lr*lr*lf*Csr);
VehParam.b11 = C1*g*lr*Csf; 
VehParam.b12 = C2*g*lf*lr*Csf; 
VehParam.a61 = etat*mu*Csf*Fzf/Js;
VehParam.a62 = etat*mu*Csf*Fzf*lf/(Js*Vx);
VehParam.a65 = -etat*mu*Csf*Fzf/Js;
VehParam.a66  = -Rs*Bs/Js;
VehParam.b16 = Rs/Js;
VehParam.Vx = Vx;
VehParam.IP = 2;


%% generate system matrices
if dim ==2 % 2 states model
    A_cont = [VehParam.a11,VehParam.a12;VehParam.a21,VehParam.a22];
    B_cont = [VehParam.b11;VehParam.b12];
    C_cont = [0 VehParam.Vx]; % single measurement
elseif dim ==4 % 4 states model
    A_cont = [VehParam.a11,VehParam.a12,0,0;VehParam.a21,VehParam.a22,0,0; 0,1,0,0; VehParam.Vx,VehParam.Vx,0,0];
    B_cont = [VehParam.b11;VehParam.b12;0;0];
    C_cont = [0 VehParam.Vx 0 0; 0 0 1 0; 0 0 0 1]; % 3 measurements
elseif dim==6 %(Comment: why are b-values in A matrix?)
    A_cont =  [VehParam.a11,VehParam.a12,0,0, VehParam.b11,0; VehParam.a21,VehParam.a22,0,0, VehParam.b12,0; 0,1,0,0,0,0; VehParam.Vx,VehParam.Vx,0,0,0,0; 0,0,0,0,0,1;VehParam.a61,VehParam.a62,0,0,VehParam.a65,VehParam.a66];
    B_cont = [0;0;0;0;0;VehParam.b16];
    C_cont = [0 VehParam.Vx 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0]; % 4 measurements
end

% convert to discrete-time model parameters
sys.F = diag(Vk(1:nrOfOutputs));
sys.E = diag(Wk(1:dim));
PlantCT = ss(A_cont,B_cont,C_cont,[]); % Continious Time Plant
PlantDT = c2d(PlantCT,options.timeStep,'impulse'); % Discretized Plant
PlantDT = ss(PlantDT.A,[PlantDT.B sys.E],PlantDT.C,0,options.timeStep); % Integrated Discretized New Plant with
A = PlantDT.A; 
B = PlantDT.B(:,1:end-size(sys.E,1));
sys.E = PlantDT.B(:,size(B,2)+1:end);
C = PlantDT.C;

% instantiate vehicle object
vehicle = linearSysDT('vehicle',A, B, zeros(length(A),1), C, options.timeStep); %initialize system


%% Parameters

params.tFinal = 3.5; %final time
params.R0 = zonotope([zeros(length(A),1),eye(length(A))]); % Initial State bounded in unity box
params.V = sys.F*zonotope([zeros(size(C,1),1),eye(size(C,1))]); % sensor noise set
params.W = sys.E*zonotope([zeros(size(sys.E,1),1),eye(size(sys.E,1))]); % disturbance set
params.u = InpSig'; %input transition


%% Simulation Settings

options.points = 1;
options.p_conf = 0.999; % probability that sample of normal distribution within specified set
options.type = 'gaussian';

% simulate result assuming Gaussian distributions
simRes = simulateRandom(vehicle, params, options);

%% obtain output values
for i=1:length(simRes.t{1})
    % create measurement noise
    v = randPoint(params.V,'gaussian',options.p_conf);
    % obtain output value
    params.y(:,i) = C*simRes.x{1}(i,:)' + v;
end

save([savepath '/' 'vehicleModel_dim_new', num2str(dim)], 'vehicle', 'params', 'options', 'simRes');


% ------------------------------ END OF CODE ------------------------------
