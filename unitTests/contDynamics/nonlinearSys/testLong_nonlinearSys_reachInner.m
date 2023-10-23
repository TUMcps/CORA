function res = testLong_nonlinearSys_reachInner
% testLong_nonlinearSys_reachInner - test if the computed
%    inner-approximation of the reachable set is correct
%
% Syntax:
%    res = testLong_nonlinearSys_reachInner
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 

% Authors:       Niklas Kochdumper
% Written:       26-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameters --------------------------------------------------------------

params.tFinal = 1;
params.R0 = interval([0.9;0],[1;0.1]);


% Reachability Settings ---------------------------------------------------

% settings for inner-approximation
options.algInner = 'scale';

options.splits = 2;
options.iter = 2;
options.orderInner = 5;
options.scaleFac = 0.95;

% settings for outer-approximation
options.timeStep = 0.001;                           
options.taylorTerms = 10;                            
options.zonotopeOrder = 50;       
options.intermediateOrder = 20;
options.errorOrder = 10;


% System Dynamics ---------------------------------------------------------

brusselator = @(x,u) [1-2*x(1) + 3/2 * x(1)^2*x(2); ...
                      x(1)-3/2*x(1)^2*x(2)];

sys = nonlinearSys(brusselator);  


% Reachability Analysis ---------------------------------------------------

% compute inner-approximation
[Rin,~] = reachInner(sys,params,options);


% Verification ------------------------------------------------------------

% Test 1: check if all points inside the computed inner-approximation
%         are located in the initial set if simulated backward in time

sysInv = nonlinearSys(@(x,u) -brusselator(x,u));

R_i = Rin.timePoint.set{end};
points = [randPoint(R_i,1000),randPoint(R_i,'all','extreme')];
points_ = zeros(size(points));

res1 = true;
for i = 1:size(points,2)

    % simulate backwards in time
    simOpts.x0 = points(:,i);
    simOpts.tFinal = params.tFinal;

    [~,x] = simulate(sysInv,simOpts);

    % check if final point is inside initial set
    p = x(end,:)';
    points_(:,i) = p;

    if ~contains(params.R0,p)
        % not a valid inner-approximation
        res1 = false;
    end
end

% figure; hold on;
% plot(params.R0);
% plot(points_(1,:),points_(2,:),'.k');

% Test 2: check if the result matches the stored one
I_saved = interval([0.665734495820338; 0.511586092476733], ...
                   [0.729556137539839; 0.565744670853688]);
I = interval(R_i);

res2 = isequal(I,I_saved,1e-3);

% combine results
res = res1 && res2;

% ------------------------------ END OF CODE ------------------------------
