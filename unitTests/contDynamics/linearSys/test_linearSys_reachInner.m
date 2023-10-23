function res = test_linearSys_reachInner()
% test_linearSys_reachInner - test if the computed inner-approximation of 
%    the reachable set is correct
%
% Syntax:
%    res = test_linearSys_reachInner()
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

% Parameter -----------------------------------------------------------

params.tFinal = 1;
params.R0 = zonotope([ones(5,1),0.1*diag(ones(5,1))]);


% Reachability Settings -----------------------------------------------

options.timeStep = 0.02;
options.zonotopeOrder = 20;


% System Dynamics -----------------------------------------------------

A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = 1;

sys = linearSys('fiveDimSys',A,B); 


% Reachability Analysis -----------------------------------------------

% compute inner-approximation
Rin = reachInner(sys,params,options);


% Verification --------------------------------------------------------

% Test 1: check if all points inside the computed inner-approximation
%         are located in the initial set if simulated backward in time

sysInv = linearSys(-A,-B); 

R_i = Rin.timePoint.set{end};
N = 20;
points = zeros(sys.dim,N);

R0 = polytope(params.R0);

res1 = true;
for i = 1:N

    % draw random point from inner-approximation
    if i < N/2
        x0 = randPoint(R_i,1,'extreme'); 
    else
        x0 = randPoint(R_i);
    end

    % simulate backwards in time
    simOpts.x0 = x0;
    simOpts.tFinal = params.tFinal;

    [~,x] = simulate(sysInv,simOpts);

    % check if final point is inside initial set
    p = x(end,:)';
    points(:,i) = p;

    if ~contains(R0,p,'exact',1e-4)
        % not a valid inner-approximation
        res1 = false;
        break
    end
end

%     figure; hold on;
%     plot(params.R0);
%     plot(points(1,:),points(2,:),'.k');


% Test 2: check if the result matches the stored one

I_saved = interval([-0.013937383819512;-0.570761541921580;0.061914997162599;-0.021873749737811;0.121801754912951], ...
                   [0.089837441984412;-0.466986716117656;0.075673885420954;-0.008114861479455;0.148868811560274]);
I = interval(R_i);

res2 = isequal(I,I_saved,1e-12);

% combine results
res = res1 && res2;
    
end
    
% ------------------------------ END OF CODE ------------------------------
