
% directory of converted spaceEx model
dirConvertedModel = [coraroot filesep 'models' filesep ...
    'SpaceExConverted' filesep 'Ranger_nonlinear_working'];

% only convert spaceEx model once
if ~exist(dirConvertedModel,'dir')
    spaceex2cora('Ranger_nonlinear_working');
end

% instantiate model
pHA = Ranger_nonlinear_working;
% correspondances of model states between spaceEx model and CORA model
% can be inferred from the generated model file
% "models\SpaceExConverted\Ranger_nonlinear_working.m"


% TODO:
% - set params (see CORA manual, Section 4): model parameters
% - set options (see CORA manual, Section 4): reachability options

% Give single initial values as a zonotope
c = [
    0; % ls; 
    0; % as; 
    0; % timer; 
    0; % turnTimer; 
    0; % arg_x; 
    0; % arg_y; 
    0; % timer; 
    0; % pos_x; 
    3; % pos_y; 
    0; % vel_x;
    0; % vel_y; 
    0; % acc_x; 
    0; % acc_y; 
    0; % yaw; 
    0; % yawVel; 
    0; % yawAcc; 
    0; % t; 
    0; % pos_x; 
    3; % pos_y; 
    0; % pos_x; 
    3]; % pos_y

G = [
    0; % ls; 
    0; % as; 
    0; % timer; 
    0; % turnTimer; 
    0; % arg_x; 
    0; % arg_y; 
    0; % timer; 
    0; % pos_x; 
    0; % pos_y; 
    0; % vel_x;
    0; % vel_y; 
    0; % acc_x; 
    0; % acc_y; 
    0; % yaw; 
    0; % yawVel; 
    0; % yawAcc; 
    0; % t; 
    0; % pos_x; 
    0; % pos_y; 
    0; % pos_x; 
    0]; % pos_y

params.tStart = 0;
params.tFinal = 10;
params.startLoc = [1;1;1;1;1];
params.R0 = zonotope(c, G);
% ...

options.guardIntersect = 'levelSet';
options.timeStep = 0.1;
options.zonotopeOrder = 50;
options.enclose = {'box'};
options.taylorTerms = 8;
% ...


% reachability analysis
R = reach(pHA,params,options);


% visualization
figure; hold on; box on;
plot(R);
