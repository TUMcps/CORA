
clear all;

% Parameters --------------------------------------------------------------

params.tFinal = 5;
R0 = interval([-0.1;0.6],[0.1;0.8]);
% original: interval([-0.5;-0.7],[0.3;0.8]);

% Reachability Settings ---------------------------------------------------

% options.timeStep = 0.01;
% options.taylorTerms = 10;
% options.zonotopeOrder = 100;
% options.intermediateOrder = 50;
% options.errorOrder = 20;
% % options.lagrangeRem.method = 'taylorModel';
% % options.lagrangeRem.simplify = 'simplify';
% 
% options.alg = 'poly';
% options.tensorOrder = 3;
% 
% options.polyZono.maxDepGenOrder = 50;
% options.polyZono.maxPolyZonoRatio = 0.01;
% options.polyZono.restructureTechnique = 'reduceFullGirard';

options.verbose = true;
options.alg = 'poly-adaptive';


% System Dynamics ---------------------------------------------------------

STB1eq = @(x,u) [-x(1)^3+x(2); -x(1)^3-x(2)^3];
sys = nonlinearSys(STB1eq);

% vector field
% [X,Y] = meshgrid(-1:0.25:1,-1:0.25:1);
% U = -X.^3 + Y;
% V = -X.^3 - Y.^3;
% figure; hold on; box on;
% quiver(X,Y,U,V);
% close;


% Reachability Analysis ---------------------------------------------------

% compute reachable set
tic;
if contains(options.alg,'poly')
    params.R0 = polyZonotope(R0);
elseif contains(options.alg,'lin')
    params.R0 = zonotope(R0);
end
R = reach(sys, params, options);
tComp = toc;
disp(['computation time (zonotope): ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOptions.points = 100;
simOptions.fracVert = 0.5;
simOptions.fracInpVert = 0.5;
simOptions.inpChanges = 1;

params.R0 = zonotope(R0);
simRes = simulateRandom(sys,params,simOptions);


% Visualization -----------------------------------------------------------
    
figure; hold on;

% plot reachable set
plot(R,[1,2],'FaceColor',[0.9290, 0.6940, 0.1250],...
    'EdgeColor',[0.9290, 0.6940, 0.1250]);

% plot initial set
plot(R0,[1,2],'w','Filled',true,'EdgeColor','k');

% plot simulation results
plot(simRes,[1,2],'b');

