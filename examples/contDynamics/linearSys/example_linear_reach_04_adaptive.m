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
%    res - boolean
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      08-October-2019
% Last update:  23-April-2020 (restructure params/options)
% Last revision:---

%------------- BEGIN CODE --------------

% System Dynamics ---------------------------------------------------------

A = [-0.7 -2; 2 -0.7];
B = 1;

sys = linearSys('sys',A,B);


% Parameter ---------------------------------------------------------------

dim = length(A);

params.tFinal = 5;
params.R0 = zonotope([[10; 5],0.5*eye(dim)]);       % initial set
params.U = zonotope([ones(dim,1),0.25*eye(dim)]);   % uncertain inputs
params.u = 0.1*[2 1 1 2 4 4 4 4 2 -2; -1 0 0 2 3 3 3 3 1 -2];

% Reachability Settings ---------------------------------------------------

options.linAlg = 'adaptive'; % use adaptive parameter tuning
options.verbose = true;

% Simulation --------------------------------------------------------------

simOpt.points = 100;
simOpt.fracVert = 0.05;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 2;

simRes = simulateRandom(sys, params, simOpt);


% Reachability Analysis ---------------------------------------------------

errs = [1,0.01];
stepssS = zeros(length(errs),1);
timesS = zeros(length(errs),1);
R = cell(length(errs),1);

% compute reachable sets for different max. allowed errors
for i=1:length(errs)
    options.error = errs(i);
    tic
    R{i} = reach(sys,params,options);
    timesS(i) = toc;
    stepssS(i) = length(R{i}.timeInterval.set);
end

% Visualization -----------------------------------------------------------

figure; hold on; box on;
projDims = [1,2];

% plot reachable set
plot(R{1},projDims,'k','EdgeColor','k');
plot(R{2},projDims,'FaceColor',[0.7,0.7,0.7],'EdgeColor',[0.7,0.7,0.7]);

% plot initial set
plot(params.R0,projDims,'w','LineWidth',1.5);

% plot simulation
plot(simRes,projDims,'b','LineWidth',0.5);

% plot unsafe set
unsafeSet = interval([2;-2],[4;2]);
plot(unsafeSet,projDims,'FaceColor',[227,114,34]/255,'Filled',true,...
    'EdgeColor','r','LineWidth',2);

% formatting
xlabel('x_1'); ylabel('x_2');
title('2D system');

res = 1;

end

%------------- END OF CODE --------------