% LECTURE 3: Reachability Analysis of Linear Systems
% required: CORA2021 - https://github.com/TUMcps/CORA
% see manual for installation of additionally required toolboxes

% Slide "Adaptive Time Steps"
% Comparison of adaptive and fixed time step sizes

close all; clear;

% System dynamics ---------------------------------------------------------

A = [-1 -4; 4 -1];
B = 1;
twoDimSys = linearSys(A,B);
dim_x = length(A);


% Model parameters --------------------------------------------------------

params.tFinal = 4;
params.R0 = zonotope([2;1],0.25*eye(dim_x));
params.U = zonotope(zeros(dim_x,1),0.025*eye(dim_x));


% Reachability analysis: fixed time step size -----------------------------

options_fixed.timeStep = 0.04;
options_fixed.taylorTerms = 5;
options_fixed.zonotopeOrder = 30;
options_fixed.linAlg = 'fromStart';

R_fixed = reach(twoDimSys, params, options_fixed);


% Reachability analysis: adaptive time step size --------------------------

options_adaptive.linAlg = 'adaptive';
options_adaptive.error = 0.015;
R_adaptive = reach(twoDimSys, params, options_adaptive);


% Visualization -----------------------------------------------------------

% colors
color_red = [227, 27, 35] ./ 255;
color_blue = [0, 92, 171] ./ 255;
% text
text_fixed = 'Fixed $\Delta t$';
text_adaptive = 'Adaptive $\Delta t$';
% line width and font size for legend
lw = 1.5;
fs = 12;

figure;

% plot: reachable sets
subplot(1,2,1); hold on; box on; axis square;
% fixed time step size
h_fixed = plot(R_fixed,[1,2],'FaceColor','w','EdgeColor',color_blue,...
    'Unify',true,'LineWidth',lw,'Filled',false);
% adaptive time step size
h_adaptive = plot(R_adaptive,[1,2],'FaceColor','w','EdgeColor',color_red,...
    'Unify',true,'LineWidth',lw,'LineStyle','--','Filled',false);
% initial set
plot(params.R0,[1,2],'k','LineWidth',lw);
% text
text(1.6,0.5,'Initial set','interpreter','latex','FontSize',fs);
text(0.75,-0.3,'Reachable sets','interpreter','latex','FontSize',fs);
% axis labels
xlabel('$x_1$','interpreter','latex','FontSize',fs);
ylabel('$x_2$','interpreter','latex','FontSize',fs);
% legend
legend([h_fixed,h_adaptive],text_fixed,text_adaptive,...
    'interpreter','latex','FontSize',fs,'Location','southeast');


% plot: time step sizes
subplot(1,2,2); hold on; box on; axis square;
% data for adaptive time step size
tVec_adaptive = query(R_adaptive,'tVec');
cumsumtVec_adaptive = [0;repelem(cumsum(tVec_adaptive(1:end-1)),2);...
    params.tFinal];
% fixed time step size
h_fixed =  plot([0,params.tFinal],...
    [options_fixed.timeStep,options_fixed.timeStep],...
    'Color',color_blue,'LineWidth',lw);
% number of steps (caution: position hardcoded)
text(3,0.032,[num2str(params.tFinal/options_fixed.timeStep) ' Steps'],...
    'interpreter','latex','FontSize',fs,'Color',color_blue);
% adaptive time step size
h_adaptive = plot(cumsumtVec_adaptive,repelem(tVec_adaptive,2),...
    'Color',color_red,'LineWidth',lw);
% number of steps (caution: position hardcoded)
text(1.5,0.08,[num2str(length(tVec_adaptive)) ' Steps'],...
    'interpreter','latex','FontSize',fs,'Color',color_red);
% axis labels
xlabel('$t$','interpreter','latex','FontSize',fs);
ylabel('$\Delta t$','interpreter','latex','FontSize',fs);
% legend
legend([h_fixed,h_adaptive],text_fixed,text_adaptive,...
    'interpreter','latex','Location','northwest','FontSize',fs);



