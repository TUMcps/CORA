function res = test_reach_poly()
% test_reach_poly - unit_test_function of nonlinear reachability analysis
%    with the conservative polynomialization approach
%
% Checks if the reachable set contains all simulated points
%
% Syntax:  
%    res = test_reach_poly()
%
% Inputs:
%    no
%
% Outputs:
%    res - boolean 

% Author:       Niklas Kochdumper
% Written:      04-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 0;

% Parameter ---------------------------------------------------------------

params.tFinal = 0.01;
params.R0 = polyZonotope(zonotope([1.4; 2.3],[0.05 0; 0 0.05]));
params.U = zonotope([1;2] + interval([-0.1;-0.1],[0.1;0.1]));


% Reachability Settings ---------------------------------------------------

options.timeStep = params.tFinal;
options.taylorTerms = 4;
options.zonotopeOrder = 10;
options.intermediateOrder = 10;
options.errorOrder = 5;

options.alg = 'poly';
options.tensorOrder = 3;


% System Dynamics ---------------------------------------------------------

f = @(x,u) [x(2)*u(1) + u(1)*u(2);
            (1-x(1))*x(2)-x(1) + u(2)^2];

sys = nonlinearSys(f);


% Reachability Analysis ---------------------------------------------------

R = reach(sys, params, options);


% Simulation --------------------------------------------------------------

% use boundary of the initial set as new initial set to get critical points 
c = center(zonotope(params.R0));
G = generators(zonotope(params.R0));

R0{1} = zonotope(c + G(:,1),G(:,2));
R0{2} = zonotope(c - G(:,1),G(:,2));
R0{3} = zonotope(c + G(:,2),G(:,1));
R0{4} = zonotope(c - G(:,2),G(:,1));

% simulation options
simOpt.points = 5000;
simOpt.fracVert = 0.1;
simOpt.fracInpVert = 0.9;
simOpt.inpChanges = 2;

% simulate the system
points = [];

for i = 1:length(R0)
   
    params.R0 = R0{i};
    simRes = simulateRandom(sys, params, simOpt);
    
    for j = 1:length(simRes.x)
       points = [points, simRes.x{j}(end,:)']; 
    end
end



% Verification ------------------------------------------------------------

% check if all points are located inside the time point reachable set
pgon = polygon(R.timePoint.set{end},12);

for i = 1:size(points,2)
   if ~in(pgon,points(:,i))
      error('Unit test failed!');
   end
end

% % visualize the set
% figure; hold on;
% plot(R.timePoint.set{end},[1,2],'r','Splits',12);
% plot(points(1,:),points(2,:),'.k');

res = 1;

%------------- END OF CODE --------------