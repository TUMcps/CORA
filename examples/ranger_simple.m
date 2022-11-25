
% directory of converted spaceEx model
dirConvertedModel = [coraroot filesep 'models' filesep ...
    'SpaceExConverted' filesep 'Ranger_simple1'];

% only convert spaceEx model once
if ~exist(dirConvertedModel,'dir')
    spaceex2cora('Ranger-simple1');
end

% instantiate model
pHA = Ranger_simple1;
% correspondances of model states between spaceEx model and CORA model
% can be inferred from the generated model file

x1 = [
    0; % pos 
    0; % vel 
    1; % acc 
    0]; % t 
x2 = [
    0; % pos 
    0; % vel 
    1; % acc 
    0]; % t

u1 = [
    1; % const_acc
    5]; % bound

u2 = [
    1; % const_acc
    5]; % bound

params.tStart = 0;
params.tFinal = 10;
params.startLoc = 1;
params.R0 = zonotope(interval(x1,x2));
params.U = zonotope(interval(u1,u2));
% ...

options.guardIntersect = 'polytope';
options.timeStep = 0.1;
options.enclose = {'box'};
options.taylorTerms = 10;
options.zonotopeOrder = 20;
options.linAlg = 'standard';
options.intersectInvariant = false;
% ...


% reachability analysis
R = reach(pHA,params,options);


% visualization
figure; hold on; box on;
plot(R);