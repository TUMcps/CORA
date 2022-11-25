% Lotka-Volterra example

clear all;

% Parameters --------------------------------------------------------------

params.tFinal = 5;
R0 = interval([0.40;0.18],[0.52;0.27]);
R0 = partition(R0,3);
R0 = R0([1,5,9]);
% original: interval([0.40;0.18],[0.52;0.27])

% Reachability Settings ---------------------------------------------------

% options.timeStep = 0.005;
% options.taylorTerms = 5;
% options.zonotopeOrder = 300;
% options.intermediateOrder = 20;
% options.errorOrder = 50;
% 
% % reachability algorithm
% options.alg = 'poly';
% options.tensorOrder = 3;
% 
% options.polyZono.maxDepGenOrder = 50;
% options.polyZono.maxPolyZonoRatio = 0.01;
% options.polyZono.restructureTechnique = 'reduceFullGirard';

options.verbose = true;
options.alg = 'poly-adaptive';


% System Dynamics ---------------------------------------------------------

f = @(x,u) [x(1)*(1.5-x(2)); -x(2)*(3-x(1))];
lv = nonlinearSys(f); 


% Reachability Analysis ---------------------------------------------------

R = cell(length(R0),1);
tic;
for i=1:length(R0)
    if contains(options.alg,'poly')
        params.R0 = polyZonotope(R0{i});
    elseif contains(options.alg,'lin')
        params.R0 = zonotope(R0{i});
    end
    R{i} = reach(lv, params, options);
end
tComp = toc;
disp(['computation time (zonotope): ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simOptions.points = 6;
simOptions.fracVert = 0.5;
simOptions.fracInpVert = 0.5;
simOptions.inpChanges = 1;

simRes = cell(length(R0),1);
for i=1:length(R0)
    params.R0 = zonotope(R0{i});
    simRes{i} = simulateRandom(lv,params,simOptions);
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;
axis([-2,15,-2,12]);

% plot reachable set
for i=1:length(R)
    plot(R{i},[1,2],'FaceColor',[0.9290, 0.6940, 0.1250],...
        'EdgeColor',[0.9290, 0.6940, 0.1250]);
end

% plot initial set
for i=1:length(R0)
    plot(R0{i},[1,2],'w','Filled',true,'EdgeColor','k');
end

% plot simulation results
for i=1:length(simRes)
    plot(simRes{i},[1,2],'b');
end

