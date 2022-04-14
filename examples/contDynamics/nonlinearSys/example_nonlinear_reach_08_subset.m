function res = example_nonlinear_reach_08_subset()
% example_nonlinear_reach_08_subset - example that demonstrates how to 
%    extract subsets of reachable sets as described in [1].
%
% Syntax:  
%    example_nonlinear_reach_08_subset()
%
% Inputs:
%    -
%
% Outputs:
%    res - 1 if example completed successfully, 0 otherwise
%
% References: 
%   [1] N. Kochdumper et al. "Utilizing Dependencies to Obtain Subsets of 
%       Reachable Sets", HSCC 2020
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyZonotope/getSubset

% Author:       Niklas Kochdumper
% Written:      28-January-2020
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% Parameter ---------------------------------------------------------------

R0 = interval([-0.2;-0.2],[0.2;0.2]) + [-1;1];
params.R0 = polyZonotope(R0);
params.tFinal = 1;


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.01;
options.taylorTerms = 4;
options.zonotopeOrder = 50;
options.intermediateOrder = 10;
options.errorOrder = 5;

options.alg = 'poly';
options.tensorOrder = 3;


% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@vanderPolEq);



% Reachability Analysis ---------------------------------------------------
     
% extract a subset of a reachable set as shown in Fig. 3 in [1]

% reachable set computation
R = reach(sys, params, options);

% extract subset from the reachable set
alpha = [0.5;0.4];                             % factor for initial point
p = getSubset(params.R0,1:2,alpha);            % initial point

Rfin = R.timePoint.set{end};
Rfin = reduce(Rfin,'girard',3);

tic
Rsubset = getSubset(Rfin,1:2,alpha);           % reachable subset
tComp = toc;

disp(['Computation time (subset extraction): ',num2str(tComp),'s']);

% simulation from the initial point (to check correctness)
paramSim.x0 = p;
paramSim.tFinal = 1;
[~,x] = simulate(sys,paramSim);

% visualization of the time interval reachable set
figure; hold on;

plot(R,[1,2],'FaceColor',[.7 .7 .7],'EdgeColor','none','Splits',4); 


% visualization of the initial and final reachable set
plot(R0,[1,2],'w','Filled',true,'EdgeColor','none');
plot(R0,[1,2],'FaceColor',[30,119,255]./255,'Filled',true,'EdgeColor','none','FaceAlpha',0.7);
plot(Rfin,[1,2],'w','Filled',true,'EdgeColor','none','Splits',10);
plot(Rfin,[1,2],'FaceColor',[30,119,255]./255,'Filled',true,'EdgeColor','none', ...
           'FaceAlpha',0.7,'Splits',10);

% visualization of the reachable subset
plot(zonotope(Rsubset),[1,2],'w','Filled',true,'EdgeColor','none');
plot(zonotope(Rsubset),[1,2],'FaceColor',[255,140,0]./255,'Filled',true,'EdgeColor','none');

% visualization of the simulation results
plot(x(end,1),x(end,2),'.k','MarkerSize',30);
plot(paramSim.x0(1),paramSim.x0(2),'Color','w','MarkerSize',20,'Marker','.');
plot(paramSim.x0(1),paramSim.x0(2),'Color',[255,140,0]./255, ...
     'MarkerSize',20,'Marker','.');

% formatting
xlim([-1.5,1.5]);
ylim([0.5 3.5]);
set(gcf,'Position',[855 434 560 375]);
box on
xlabel('$x_1$','FontSize',20,'interpreter','latex');
ylabel('$x_2$','FontSize',20,'interpreter','latex');



% Falsification -----------------------------------------------------------

% compute a falsifying trajectory as described in Sec. 4.1 in [1]


% constraint defining the specification
hs = halfspace(-[1,2],-6.4);

% construct objective function
Rproj = -hs.c' * Rfin;

b = -hs.d - Rproj.c - sum(abs(Rproj.Grest));
x = sym('x',[size(Rproj.expMat,1),1]);

objSym = -b;

for i = 1:length(Rproj.G)
   temp = x(1)^Rproj.expMat(1,i);
   for j = 2:length(x)
       temp = temp*x(j)^Rproj.expMat(j,i);
   end
   objSym = objSym + Rproj.G(i)*temp;
end

objFun = matlabFunction(-objSym,'Vars',{x});

% compute most critical initial point by optimization (see (28) in [1])
lb = -ones(size(Rproj.expMat,1),1);
ub = ones(size(Rproj.expMat,1),1);

x0 = zeros(size(Rproj.expMat,1),1);

tic
alpha = fmincon(objFun,x0,[],[],[],[],lb,ub);
p = getSubset(params.R0,1:2,alpha);

% compute the falsifying trajectory
paramSim.x0 = p;
[~,xTraj] = simulate(sys,paramSim);
tComp = toc;

disp(['Computation time (falsification): ',num2str(tComp),'s']);

figure
hold on

% visualize constraint
xlim([-1.5,1.5]);
ylim([0.5 3.5]);
plot(hs,[1,2],'r','FaceAlpha',0.5,'EdgeColor','none');

% visualize time interval reachable set
plot(R,[1,2],'FaceColor',[.7 .7 .7],'EdgeColor','none','Splits',4); 

% visualize initial and final reachable set
plot(params.R0,[1,2],'w','Filled',true,'EdgeColor','none');
plot(Rfin,[1,2],'w','Filled',true,'EdgeColor','none','Splits',10);
plot(params.R0,[1,2],'FaceColor',[255,140,0]./255,'Filled',true,'EdgeColor', ...
                    'none','FaceAlpha',0.7);
plot(Rfin,[1,2],'FaceColor',[30,119,255]./255,'Filled',true,'EdgeColor','none', ...
                    'FaceAlpha',0.7,'Splits',10);

% visualize falsifying trajectory
plot(xTraj(:,1),xTraj(:,2),'k');
plot(xTraj(1,1),xTraj(1,2),'.k','MarkerSize',15);
plot(xTraj(end,1),xTraj(end,2),'.k','MarkerSize',15);

% formatting
set(gcf,'Position',[855 434 560 375]);
box on
xlabel('$x_1$','FontSize',20,'interpreter','latex');
ylabel('$x_2$','FontSize',20,'interpreter','latex');



% Optimization ------------------------------------------------------------

% optimization over reachable sets as described in Sec. 4.2 in [1]


% objective function: maximize volume of initial set
fun = @(x) -(x(2)-x(1))*(x(4)-x(3));

A = [1 -1 0 0;0 0 1 -1];
b = [0;0];

% solve optimization problem (see (24) in [1])
tic
x = fmincon(fun,[-1;1;-1;1],A,b,[],[],-ones(4,1),ones(4,1),@conFun);
topt1 = toc;
disp(['Computation time (optimization): ',num2str(topt1)]);

% compute optimized final reachable set
int = interval([x(1);x(3)],[x(2);x(4)]);

Rfin_ = getSubset(Rfin,1:2,int);

% compute optimized initial set
R0_ = zonotope(getSubset(params.R0,1:2,int));

% visualize constraint
figure 
hold on

xlim([-1.5,1.5]);
ylim([0.5 3.5]);
plot(hs,[1,2],'r','FaceAlpha',0.5,'EdgeColor','none');

% visualize time interval reachable set
plot(R,[1,2],'FaceColor',[.7 .7 .7],'EdgeColor','none','Splits',4); 

% visualize initial and final reachable set
plot(R0,[1,2],'w','Filled',true,'EdgeColor','none');
plot(R0,[1,2],'FaceColor',[30,119,255]./255,'Filled',true,'EdgeColor','none','FaceAlpha',0.7);
plot(Rfin,[1,2],'w','Filled',true,'EdgeColor','none','Splits',10);
plot(Rfin,[1,2],'FaceColor',[30,119,255]./255,'Filled',true,'EdgeColor','none', ...
           'FaceAlpha',0.7,'Splits',10);

% visualize optimized initial and final set
plot(Rfin_,[1,2],'w','Filled',true,'EdgeColor','none');
plot(Rfin_,[1,2],'FaceColor',[255,140,0]./255,'Filled',true,'EdgeColor','none','FaceAlpha',0.7);
plot(R0_,[1,2],'w','Filled',true,'EdgeColor','none');
plot(R0_,[1,2],'FaceColor',[255,140,0]./255,'Filled',true,'EdgeColor','none','FaceAlpha',0.7);


% formatting
set(gcf,'Position',[855 434 560 375]);
box on
xlabel('$x_1$','FontSize',20,'interpreter','latex');
ylabel('$x_2$','FontSize',20,'interpreter','latex');

res = 1;

end


% Auxiliary Functions -----------------------------------------------------

function [c,ceq] = conFun(x)
% constrained function for the optimization problem (see (31) in [1])

    % compute maximum
    a = 0.05 - 0.2 * x(4);
    b = -0.22;
    
    alpha_1 = -a/(2*b);
    alpha_2 = x(4);
    
    % compute value for the constraint
    c = 0.05*alpha_1 + 0.66*alpha_2 - 0.22*alpha_1^2 - ...
        0.2*alpha_1*alpha_2 + 0.04;
    
    ceq = [];
end
    
%------------- END OF CODE -------------