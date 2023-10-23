function res = example_hybrid_reach_04_spacecraft()
% example_hybrid_reach_04_spacecraft - spacecraft rendevous example as
%                                      described in [1] solved with the 
%                                      method from [2]
%
% Syntax:
%    res = example_hybrid_reach_04_spacecraft()
%
% Inputs:
%    no
%
% Outputs:
%    res - 1 if example completed succesfully, 0 otherwise
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission (Benchmark proposal)"  
%   [2] N. Kochdumper et al. "Reachability Analysis for Hybrid Systems with 
%       Nonlinear Guard Sets", HSCC 2020

% Authors:       Niklas Kochdumper
% Written:       23-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Parameters --------------------------------------------------------------

% problem description
R0 = zonotope([[-900; -400; 0; 0],diag([25;25;0;0])]);

params.R0 = R0;                                    % initial set
params.tFinal = 200;                               % final time
params.startLoc = 1;                               % initial location


% Reachability Options ----------------------------------------------------

% algorithm options
options.timeStep{1} = 2e-1;
options.timeStep{2} = 5e-2;

options.zonotopeOrder = 40; 
options.taylorTerms = 3;

options.intermediateOrder = 2;
options.errorOrder = 5;

options.tensorOrder = 2;
options.alg = 'lin';

% guard intersection method
options.guardIntersect = 'levelSet';


% System Dynamics ---------------------------------------------------------

HA = spacecraft_levelSet();


% Reachability Analysis ---------------------------------------------------

tic
R = reach(HA,params,options);
tComp = toc;

disp(['computation time: ',num2str(tComp)]);


% Simulation --------------------------------------------------------------

simRes = simulateRandom(HA,params); 


% Visualization -----------------------------------------------------------

% Plot 1: position space

figure; hold on; grid on; box on

% plot line-of-sight cone
phi = [0:0.01:2*pi,0];
x = 100*cos(phi);
y = 100*sin(phi);
polyCirc = polyshape([x;y]');
poly = polyshape([-100,-60;0,0;-100,60]);
poly = intersect(poly,polyCirc);
plot(poly,'FaceColor',colorblind('r'),'EdgeColor','none','FaceAlpha',0.5);

% plot reachable set
plot(R);

% plot initial set
plot(params.R0,[1,2],'k','FaceColor','w');

% plot simulation
plot(simRes,[1,2]);

% plot guard set
plot(x,y,'k');

% formatting
xlabel('$x~[m]$','interpreter','latex','FontSize',15);
ylabel('$y~[m]$','interpreter','latex','FontSize',15);
axis equal
xlim([-1000,200]);
ylim([-450,150]);


% Plot 2: velocity space

figure ; hold on; grid on; box on

% plot velocity constraint
phi = [0:0.01:2*pi,0];
x = 3.3*cos(phi);
y = 3.3*sin(phi);
poly = polyshape([x;y]');
plot(poly,'FaceColor',colorblind('r'),'EdgeColor','none','FaceAlpha',0.5);

% plot reachable set for location 2
R2 = find(R,'location',2);
plot(R2,[3,4]);

% plot simulation for location 2
for i = 1:length(simRes)
    for j=1:length(simRes(i).loc)
        if simRes(i).loc(j) == 2
            x = simRes(i).x{j};
            plot(x(:,3),x(:,4),'Color',colorblind('y'));
        end
    end
end

% formatting
xlabel('$v_x~[m/min]$','interpreter','latex','FontSize',15);
ylabel('$v_y~[m/min]$','interpreter','latex','FontSize',15);
axis equal
xlim([-4,4]);
ylim([-4,4]);

res = true;

% ------------------------------ END OF CODE ------------------------------
