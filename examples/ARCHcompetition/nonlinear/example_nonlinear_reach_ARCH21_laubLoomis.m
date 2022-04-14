function example_nonlinear_reach_ARCH21_laubLoomis()
% example_nonlinear_reach_ARCH21_laubLoomis - example of 
% nonlinear reachability analysis
%
% Syntax:  
%    example_nonlinear_reach_ARCH21_laubLoomis
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% 
% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      27-March-2019
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------


dim=7;

% Parameter ---------------------------------------------------------------

params.tFinal = 20;
x0 = [1.2; 1.05; 1.5; 2.4; 1; 0.1; 0.45];


% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 200;
options.intermediateOrder = 20;
options.taylorTerms = 4;

options.lagrangeRem.simplify = 'optimize';



% Reachability Analysis ---------------------------------------------------

R = cell(3,1);

% Initial set W = 0.1
W = 0.1;
params.R0 = polyZonotope(x0,W*eye(dim));

options1 = options;
options1.timeStep = 0.02;
options1.errorOrder = 5;
options1.alg = 'poly';
options1.tensorOrder = 3;

polyZono.maxPolyZonoRatio = 0.1;
polyZono.restructureTechnique = 'reduceFullGirard';
polyZono.maxDepGenOrder = 8;
options1.polyZono = polyZono;

sys = nonlinearSys('sevenDim1',@laubLoomis); 

tic
R{1} = reach(sys, params, options1);
tComp = toc;

width = 2*rad(interval(project(R{1}.timePoint.set{end},4)));
disp(['computation time of reachable set (W=0.1): ',num2str(tComp)]);
disp(['width of final reachable set (W=0.1): ',num2str(width)]);
disp(' ');


% Initial set W = 0.05
W = 0.05;
params.R0 = zonotope([x0,W*eye(dim)]);

options2 = options;
options2.timeStep = 0.025;
options2.errorOrder = 1;
options2.alg = 'lin';
options2.tensorOrder = 3;

sys = nonlinearSys('sevenDim2',@laubLoomis); 

tic
R{2} = reach(sys, params, options2);
tComp = toc;

width = 2*rad(interval(project(R{2}.timePoint.set{end},4)));
disp(['computation time of reachable set (W=0.05): ',num2str(tComp)]);
disp(['width of final reachable set (W=0.05): ',num2str(width)]);
disp(' ');


% Initial set W = 0.01
W = 0.01;
params.R0 = zonotope([x0,W*eye(dim)]);

options3 = options;
options3.timeStep = 0.1;
options3.errorOrder = 1;
options3.tensorOrder = 2;
options3.alg = 'lin';

sys = nonlinearSys('sevenDim3',@laubLoomis); 

tic
R{3} = reach(sys, params, options3);
tComp = toc;

width = 2*rad(interval(project(R{3}.timePoint.set{end},4)));
disp(['computation time of reachable set (W=0.01): ',num2str(tComp)]);
disp(['width of final reachable set (W=0.01): ',num2str(width)]);
disp(' ');


% Visualization -----------------------------------------------------------

colors = {'b',[0,0.6,0],'r'};

figure; hold on; box on;

% plot results over time
for k = 1:length(R)
    plotOverTime(R{k},4,'FaceColor',colors{k},'EdgeColor','none');
end

xlabel('t')
ylabel('x_4')
axis([0,20,1.5,5])


%------------- END OF CODE --------------