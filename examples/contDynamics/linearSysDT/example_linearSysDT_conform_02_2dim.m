function completed = example_linearSysDT_conform_02_2dim
% example_linearSysDT_conform_02_2dim - example for conformance synthesis
%    of linear discrete-time systems from a unit test; shows the
%    solution of the linearSysDT class for a two-dimensional example
%    from Sec. 7.1 of [1].
%
% Syntax:
%    completed = example_linearSysDT_conform_02_2dim
%
% Inputs:
%    -
%
% Outputs:
%    completed - true/false
%
% Reference:
%    [1] T. Alamo, J. M. Bravo, and E. F. Camacho. Guaranteed
%        state estimation by zonotopes. Automatica, 41(6):1035–1043, 2005.

% Authors:       Matthias Althoff
% Written:       26-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


%% Parameters 

params.tFinal = 5; %final time
params.R0 = zonotope(zeros(2,1),eye(2)); % initial set template for conformance
params.V = zonotope(0,1); % sensor noise set
params.W = [-6; 1]*zonotope([0,1]); % disturbance set
y = [0.79; 5.00; 4.35; 1.86; -0.11; -1.13]'; % measurement vector
delta = [0.1; 0.2; 0.1; 0; -0.1; -0.2]'; % deviation of measurement vector
params.testSuite = [trajectory(zeros(2,6), [1;1], y, [], 1); ... % trajectory 1
    trajectory(zeros(2,6), [1;1], y + delta, [], 1); ... % trajectory 2
    trajectory(zeros(2,6), [1;1], y - delta, [], 1)]; % trajectory 3

%% Algorithmic Settings 
options = struct;


%% System Dynamics 

A = [0 -0.5; 1 1];
B = [1; 1];
c = zeros(2,1);
C = [-2 1];
E = eye(dim(params.W));
F = eye(dim(params.V));
twoDimSys = linearSysDT('twoDimSys',A,B,c,C,[],[],E,F,1); 


%% Reachset Synthesis
[params_conform, results] = conform(twoDimSys,params,options);
union_y_a = results.unifiedOutputs;

%% Compute reachable set using obtained parameters
options.zonotopeOrder = inf;
R = reach(twoDimSys, params_conform, options);

%% Plot trajectories and reachable sets

figure; hold on;
% plot time elapse
plotOverTime(R,1,'LineWidth',2);
% plot unified trajectories
for k = 1:size(union_y_a,2)
    plot(k-1, squeeze(union_y_a(1,k,:)),'k','Marker','.','LineStyle', 'none');
end

% label plot
xlabel('t');
ylabel('y');

completed = true;

% ------------------------------ END OF CODE ------------------------------
