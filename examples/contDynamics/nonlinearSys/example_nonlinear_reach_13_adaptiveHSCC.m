function res = example_nonlinear_reach_13_adaptiveHSCC
% example_nonlinear_reach_13_adaptiveHSCC - example for nonlinear
%    reachability analysis using adaptive parameter tuning,
%    reproducing results from [1]
%
% Syntax:
%    res = example_nonlinear_reach_13_adaptiveHSCC
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References:
%    [1] M. Wetzlinger, A. Kulmburg, M. Althoff. "Adaptive Parameter Tuning
%        for Reachability Analysis of Nonlinear Systems", HSCC 2021.

% Authors:       Mark Wetzlinger
% Written:       02-February-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% plotting
makePlots = true;
fontsize = 17;

tab2 = zeros(2,4);

% PRDE20 ------------------------------------------------------------------

% init system
dim_x = 3;
params.tFinal = 100;
params.R0 = zonotope(interval([9.5;0.01;0.01],[10;0.01;0.01]));
params.U = zonotope(0);
sys = nonlinearSys(@prodDes,dim_x,1);

% linearization approach
options.alg = 'lin-adaptive';
options.verbose = true;

adapTime = tic;
[R,~,opt] = reach(sys,params,options);
endset = R.timePoint.set{end};
tComp = toc(adapTime);
tab2(1,1) = tComp;
tab2(1,2) = volume(interval(endset));

% simulation
simOpt.points = 10;           % number of initial points
simOpt.fracVert = 0.2;        % fraction of vertices initial set

simRes = simulateRandom(sys,params,simOpt);

% figures
if makePlots
% 1. reachable sets and simulation
fig6 = figure; rows = 1; cols = 2;
projDims = {[1,2],[2,3]};
sysaxis = {[-0.5,10.5,-0.5,6.5],[-0.5,6.5,-0.5,10.5]};
sysxticks = {[0,5,10],[0,3,6]};
sysyticks = {[0,3,6],[0,5,10]};
for p=1:ceil(sys.dim/2)
    subplot(rows,cols,p); hold on; box on;
    useCORAcolors("CORA:contDynamics")
    plot(R,projDims{p});
    plot(R(1).R0);
    plot(simRes,projDims{p});
    if p == 1; plot(params.R0,projDims{p},'r','LineWidth',1.5);
    else;      scatter(0.01,0.01,3,'r','filled');
    end

    axis(sysaxis{p});
    ax = gca; ax.FontSize = 13;
    xlabel(['$x_',num2str(projDims{p}(1)),'$'],...
        'Fontsize',fontsize,'interpreter','latex');
    ylabel(['$x_',num2str(projDims{p}(2)),'$'],...
        'Fontsize',fontsize,'interpreter','latex');
    xticks(sysxticks{p});
    yticks(sysyticks{p});
end
sgtitle("Figure 6");

% time step size
fig7a = figure; hold on; box on;
title("Figure 7a");
useCORAcolors("CORA:default")
tVec = query(R,'tVec');
cumsumtVec = cumsum(tVec);
tVecSteps = [0;repelem(cumsumtVec(1:end-1),2);cumsumtVec(end)];
plot(tVecSteps,repelem(tVec,2));
% axes and labels
ax = gca; ax.FontSize = 11;
xlabel('t','FontSize',fontsize,'interpreter','latex');
ylabel('$\Delta t$','FontSize',fontsize,'interpreter','latex');        

% taylor terms (Rlin and Rerr)
fig7b = figure; hold on; box on;
title("Figure 7b");
useCORAcolors("CORA:default")
plot(tVecSteps,repelem(opt.tt_lin,2));
plot(tVecSteps,repelem(opt.tt_err,2));
axis([0,params.tFinal,0,max([opt.tt_lin;opt.tt_err])+1]);
ax = gca; ax.FontSize = 11;
legend('$\eta_{lin}$','$\eta_{abs}$','Location','northeast','Orientation','horizontal',...
    'FontSize',fontsize-2,'interpreter','latex');
legend box off;
xlabel('t','FontSize',fontsize,'interpreter','latex');
ylabel('$\eta$','FontSize',fontsize,'interpreter','latex');

% zonotope order
fig7c = figure; hold on; box on;
title("Figure 7c");
useCORAcolors("CORA:default")
fullzonorderRtp = sum(opt.zonordersRtp,2);
plot(tVecSteps,repelem(fullzonorderRtp,2));
% axes and labels
axis([0,params.tFinal,0,ceil(1.1*max(fullzonorderRtp))]);
ax = gca; ax.FontSize = 11;
xlabel('t','FontSize',fontsize,'interpreter','latex');
if strcmp(options.alg,'lin')
    ylabel('$\rho$','FontSize',fontsize,'interpreter','latex');
elseif strcmp(options.alg,'poly')
    ylabel('$\rho$','FontSize',fontsize,'interpreter','latex');
end
end

% Poly. approach
options.alg = 'poly-adaptive';

adapTime = tic;
[R,~,~] = reach(sys,params,options);
endset = R.timePoint.set{end};
tComp = toc(adapTime);
tab2(2,1) = tComp;
tab2(2,2) = volume(interval(endset));


% LALO20 ------------------------------------------------------------------

% init system
dim_x = 7;
params.tFinal = 20;
x0 = [1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45];
W = 0.05; % other options: 0.10, 0.05, 0.01
params.R0 = zonotope([x0,W*eye(dim_x)]);
params.U = zonotope(0);
sys = nonlinearSys(@laubLoomis,dim_x,1);

% Lin. approach
options.alg = 'lin-adaptive';
% options.verbose = true;

adapTime = tic;
[R,~,~] = reach(sys,params,options);
endset = R.timePoint.set{end};
tComp = toc(adapTime);
tab2(1,3) = tComp;
tab2(1,4) = max(2*rad(interval(endset)));

% Poly. approach
options.alg = 'poly-adaptive';

adapTime = tic;
[R,~,opt] = reach(sys,params,options);
endset = R.timePoint.set{end};
tComp = toc(adapTime);
tab2(2,3) = tComp;
tab2(2,4) = max(2*rad(interval(endset)));

if makePlots
% time step size
fig8 = figure; sgtitle("Figure 8");
subplot(2,1,1); hold on; box on;
useCORAcolors("CORA:default")
tVec = query(R,'tVec');
cumsumtVec = cumsum(tVec);
tVecSteps = [0;repelem(cumsumtVec(1:end-1),2);cumsumtVec(end)];
plot(tVecSteps,repelem(tVec,2));
% axes and labels
% axis([0,params.tFinal,0.9*min(tVec),1.1*max(tVec)]);
ax = gca; ax.FontSize = 11;
xlabel('t','FontSize',fontsize,'interpreter','latex');
ylabel('$\Delta t$','FontSize',fontsize,'interpreter','latex');


% zonotope order
subplot(2,1,2); hold on; box on;
useCORAcolors("CORA:default")
fullzonorderRtp = sum(opt.zonordersRtp,2);
plot(tVecSteps,repelem(fullzonorderRtp,2));
% axes and labels
axis([0,params.tFinal,0,ceil(1.1*max(fullzonorderRtp))]);
ax = gca; ax.FontSize = 11;
xlabel('t','FontSize',fontsize,'interpreter','latex');
if strcmp(options.alg,'lin')
    ylabel('$\rho$','FontSize',fontsize,'interpreter','latex');
elseif strcmp(options.alg,'poly')
    ylabel('$\rho$','FontSize',fontsize,'interpreter','latex');
end
end

end

% ------------------------------ END OF CODE ------------------------------
