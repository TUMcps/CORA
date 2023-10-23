function res = test_nonlinearSysDT_reach_04_uTransVec
% test_nonlinearSysDT_reach_04_uTransVec - unit test for nonlinear 
%    discrete time reachability analysis from [1, Sec.6] using uTransVec
%
% Syntax:
%    res = test_nonlinearSysDT_reach_04_uTransVec
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Reference:
%    [1] J.M. Bravo, Robust MPC of constrained discrete-time
%        nonlinear systems based on approximated reachable sets, 2006.

% Authors:       Laura Luetzow
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 0.075;
params.R0 = zonotope([-0.15;-45],diag([0.005;3]));
params.U = zonotope(zeros(2,1),diag([0.1;2]));
params.uTransVec = [0.1*ones(2,2) 0.1*rand(2, 1) 0.1*ones(2,3)]; 

% Reachability Settings ---------------------------------------------------

options.zonotopeOrder = 100;
options.errorOrder = 5;
options.tensorOrder = 3;


% System Dynamics ---------------------------------------------------------

% sampling time
dt = 0.015;

fun = @(x,u) cstrDiscr(x,u,dt);

sysDisc = nonlinearSysDT('stirredTankReactor',fun,dt);


% Reachability Analysis ---------------------------------------------------

R = reach(sysDisc,params,options);


% Simulation and Verification ---------------------------------------------

nr_points = 200;
for s = 1:nr_points
    params.x0 = randPoint(params.R0,1,'extreme');
    params.u = params.uTransVec + randPoint(params.U, 6,'extreme');
    [~,x,~,~] = simulate(sysDisc,params);
    for k=1:size(x,1)
        if ~contains(R.timePoint.set{k}, x(k,:)')
            res = false;
            return
        end
    end
end

res = true;
return 

% ------------------------------ END OF CODE ------------------------------
