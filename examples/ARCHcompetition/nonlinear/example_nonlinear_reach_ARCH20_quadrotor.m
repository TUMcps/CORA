function example_nonlinear_reach_ARCH20_quadrotor()
% example_nonlinear_reach_ARCH20_quadrotor - quadrotor example from the
%                                            ARCH20 competition
%
% Syntax:  
%    example_nonlinear_reach_ARCH20_quadrotor
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean  
% 
% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      15-April-2017
% Last update:  30-May-2018
%               13-March-2019
% Last revision:---


%------------- BEGIN CODE --------------

dim = 12;

% Parameter ---------------------------------------------------------------

params.tFinal = 5; 
params.U = zonotope([1;0;0]);

x0 = zeros(dim,1);


% Reachability Settings ---------------------------------------------------

options.timeStep = 0.1;
options.taylorTerms = 4;
options.zonotopeOrder = 50;
options.intermediateOrder = 5;
options.errorOrder = 1;

options.alg = 'lin';
options.tensorOrder = 2;
 

% System Dynamics ---------------------------------------------------------

sys = nonlinearSys(@quadrocopterControlledEq);


% Reachability Analysis ---------------------------------------------------

R = cell(3,1);

% case dx = 0.8
params.R0 = zonotope(x0,0.8*diag([ones(6,1);zeros(6,1)]));

tic
R{1} = reach(sys, params, options);
tComp1 = toc;


% case dx = 0.4
params.R0 = zonotope(x0,0.4*diag([ones(6,1);zeros(6,1)]));

tic
R{2} = reach(sys, params, options);
tComp2 = toc;


% case dx = 0.1
params.R0 = zonotope(x0,0.1*diag([ones(6,1);zeros(6,1)]));

tic
R{3} = reach(sys, params, options);
tComp3 = toc;


% Verification ------------------------------------------------------------

goal = interval(0.98,1.02);
t = options.timeStep;
res = 1;

for i = 1:length(R{2}.timeInterval.set)
    
    int = interval(project(R{2}.timeInterval.set{i},3));

    % check if height is below 1.4 for all times
    if supremum(int) >= 1.4
       res = 0;
       break;
    end

    % check if height is above 0.9 after 1 second
    if t > 1 && infimum(int) < 0.9
       res = 0;
       break;
    end
end

% check if final reachable set is in goal region
if ~in(goal,project(R{2}.timePoint.set{end},3))
   res = 0; 
end

disp(' ');
disp(['verified: ',num2str(res)]);

% case dx = 0.8
disp(['computation time (dx = 0.8): ',num2str(tComp1)]);

% case dx = 0.4
disp(['computation time (dx = 0.4): ',num2str(tComp2)]);

% case dx = 0.1
disp(['computation time (dx = 0.1): ',num2str(tComp3)]);



% Visualization -----------------------------------------------------------

figure; hold on; box on;
colors = {'b',[0,0.6,0],'r'};

% plot results over time
for k = 1:length(R)
    plotOverTime(R{k},3,'FaceColor',colors{k},'EdgeColor','none');
end
    
% label plot
xlabel('t');
ylabel('x_3');
axis([0,5,-0.8,1.5]);

%------------- END OF CODE --------------