function res = example_nonlinear_reach_09_drum()
% example_nonlinear_reach_09_drum - example for nonlinear reachability
%    analysis taken from [1]
%
% Syntax:
%    res = example_nonlinear_reach_09_drum
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
% 
% References:
%   [1] El-Guindy, Ahmed, Dongkun Han, and Matthias Althoff.
%       "Formal analysis of drum-boiler units to maximize
%       the load-following capabilities of power plants."
%       IEEE Transactions on Power Systems 31.6 (2016): 4691-4702.

% Authors:       Ahmed El-Guindy, Matthias Althoff, Mark Wetzlinger
% Written:       22-May-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Parameter ---------------------------------------------------------------

params.tStart = 0;   % Start time
params.tFinal = 300; % Final time

% Initial set for a change of 40MW from 110MW towards 70MW
x0 = [5.5e5; 19.06; 0.01027; 3.687; 967.7; 9.241e5; ...
        22e6; 39.88; 39.88; 14.44; 6.382; 6.382]; 

%Initial set for a change of 40MW from 70MW towards 110MW
% x0 = [5.5e5; 18.06; 0.01427; 2.687; 867.7; 8.241e5; ...
%         18e6; 39.88; 39.88; 20.44; 6.382; 6.382];

% Mapping of the reachable set to obtain the water level of the drum
C = [1e-5       0          0           0         0 0 0       0 0 0 0 0;
     0          0          0           0         0 0 4.75e-6 0 0 0 0 0;
     0          0          0           0         0 0 0       0 0 0 1 0;
    -4.7523e-04 6.8027e+01 1.9778e+04 6.8027e+01 0 0 0       0 0 0 0 0]; 

% Uncertain initial sets
Bound = diag([10000,0.5,0.001,0.5,0,0,100000,2,0,1,1,0]);

% Initial sets for reachability analysis
params.R0 = zonotope([x0,Bound]); 
R0Level   = C*params.R0 + [ 0;0;0 ;(-1.4778e3)];

% Initial set for the input variables
uTrans = [70; 0; 0];
Bound_u = diag([0,0.01,5]);
params.U = zonotope([uTrans,Bound_u]);


% Reachability Settings ---------------------------------------------------

options.timeStep      = 1;   % Time step size for reachable set computation
options.taylorTerms   = 20;
options.zonotopeOrder = 500; % Zonotope order

% Linearization options
options.tensorOrder   = 2;
options.alg           = 'lin';
options.maxError      = 1.0e+02 * ones(12,1);


% System Dynamics ---------------------------------------------------------

drumSys = nonlinearSys(@drum);


% Reachability Analysis ---------------------------------------------------

% Compute reachable set 
Rcont = reach(drumSys, params, options);

% Mapping to get the rechable set of the water level 
R = query(Rcont,'reachSet');
RcontY = cell(length(R),1);
for i=1:length(R)
    RcontY{i} = C*R{i} + [0; 0; 0; -1.4778e3];
end


% Simulation --------------------------------------------------------------

%random simulations
simOpt.points = 60;
simRes = simulateRandom(drumSys, params, simOpt);

% obtain level from states
level = cell(length(simRes),1);
for iRun = 1:length(simRes)
    level{iRun} = -4.7523e-04*simRes(iRun).x{1}(:,1) + ...
        6.8027e+01*simRes(iRun).x{1}(:,2) + ...
        6.8027e+01*simRes(iRun).x{1}(:,4) + ...
        1.9778e+04*simRes(iRun).x{1}(:,3) - 1.4778e3;
end 


% Visualization -----------------------------------------------------------

figure;   

% Power demand vs level
subplot(2,2,1); hold on;
for i=1:length(RcontY)    
    ZprojY = project(RcontY{i},[2 4]);
    ZprojY = reduce(ZprojY,'girard',10);
    plot(ZprojY,[1 2],'FaceColor',CORAcolor('CORA:reachSet'));
end
plot(R0Level,[2 4],'k','FaceColor',CORAcolor("CORA:initialSet"));
for i=1:length(simRes)
    plot(simRes(i).x{1}(:,7)*4.75e-6,level{i},'Color',CORAcolor("CORA:simulations"));
end

% Pressure vs level
subplot(2,2,2); hold on;
for i=1:length(RcontY)    
    ZprojY = project(RcontY{i},[1 4]);
    ZprojY = reduce(ZprojY,'girard',10);
    plot(ZprojY,[1 2],'FaceColor',CORAcolor('CORA:reachSet'));
end
plot(R0Level,[1 4],'k','FaceColor',CORAcolor("CORA:initialSet"));
for i=1:length(simRes)
    plot(simRes(i).x{1}(:,1)*1e-5,level{i},'Color',CORAcolor("CORA:simulations"));
end

% Power demand vs Pressure
subplot(2,2,3); hold on;
for i=1:length(RcontY)    
    ZprojY = project(RcontY{i},[2 1]);
    ZprojY = reduce(ZprojY,'girard',10);
    plot(ZprojY,[1 2],'FaceColor',CORAcolor('CORA:reachSet'));
end
plot(R0Level,[2 1],'k','FaceColor',CORAcolor("CORA:initialSet"));
for i=1:length(simRes)
    plot(simRes(i).x{1}(:,7)*4.75e-6,simRes(i).x{1}(:,1)*1e-5,...
       'Color',CORAcolor("CORA:simulations"));
end

% Steam quality vs Steam volume
subplot(2,2,4); hold on;
projDims=[3 4];
plot(Rcont,projDims,'FaceColor',CORAcolor('CORA:reachSet'));
plot(params.R0,projDims,'k','FaceColor',CORAcolor("CORA:initialSet"));
for i=1:length(simRes)
    plot(simRes(i).x{1}(:,projDims(1)),simRes(i).x{1}(:,projDims(2)),...
        'Color',CORAcolor("CORA:simulations"));
end
xlabel("x_" + projDims(1)); ylabel("x_" + projDims(2));

    
% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
