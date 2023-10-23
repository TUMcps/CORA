function res = example_linear_reach_04_adaptive
% example_linear_reach_04_adaptive - example for adaptive parameter tuning
%                                    for linear time-invariant systems
%
% Syntax:
%    res = example_linear_reach_04_adaptive
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       08-October-2019
% Last update:   23-April-2020 (restructure params/options)
%                07-November-2022 (change colors)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% System Dynamics ---------------------------------------------------------

A = [-0.7 -2; 2 -0.7];
B = 1;
sys = linearSys('sys',A,B);


% Parameters --------------------------------------------------------------

dim = length(A);

params.tFinal = 5;
params.R0 = zonotope([[10; 5],0.5*eye(dim)]);       % initial set
params.U = zonotope([ones(dim,1),0.25*eye(dim)]);   % uncertain inputs
params.u = 0.1*[2 1 1 2 4 4 4 4 2 -2; -1 0 0 2 3 3 3 3 1 -2];


% Reachability Settings ---------------------------------------------------

options.linAlg = 'adaptive'; % use adaptive parameter tuning


% Simulation --------------------------------------------------------------

simOpt.points = 100;
simOpt.fracVert = 0.05;
simRes = simulateRandom(sys, params, simOpt);


% Reachability Analysis ---------------------------------------------------

errs = [1,0.01];
stepssS = zeros(length(errs),1);
timesS = zeros(length(errs),1);
R = cell(length(errs),1);

% compute reachable sets for different max. allowed errors
for i=1:length(errs)
    options.error = errs(i);
    tic;
    R{i} = reach(sys,params,options);
    timesS(i) = toc;
    stepssS(i) = length(R{i}.timeInterval.set);
end

% Visualization -----------------------------------------------------------

figure; hold on; box on; legend();
projDims = [1,2];

% plot unsafe set
unsafeSet = specification(interval([2;-2],[4;2]));
plot(unsafeSet,projDims,'DisplayName', 'Unsafe Set');

% plot reachable set
useCORAcolors("CORA:contDynamics", 2)
plot(R{1}, projDims, 'DisplayName', sprintf("Reachable set: error=%.2f", errs(1))); 
plot(R{2}, projDims, 'DisplayName', sprintf("Reachable set: error=%.2f", errs(2))); 

% plot initial set
plot(R{1}.R0,projDims, 'DisplayName', 'Initial set');

% plot simulation
plot(simRes,projDims, 'DisplayName', 'Simulations');

% formatting
xlabel('x_1'); ylabel('x_2');

res = true;

end

% ------------------------------ END OF CODE ------------------------------
