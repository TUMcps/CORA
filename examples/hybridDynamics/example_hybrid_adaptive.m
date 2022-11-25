% FILE AUTOMATICALLY GENERATED ON, 16-Apr-2021 10:08:30
function [] = example_hybrid_adaptive()


% System Dynamics --------------------------------------------------------- 

hybridSys = bball();


% Parameter ---------------------------------------------------------------

params.tFinal = 1.7;
params.R0 = zonotope([1;0], diag([0.05,0.05]));
params.startLoc = 1;
params.finalLoc = 0;
params.U{1} = zonotope(interval(0, 0));


% Reachability Analysis ---------------------------------------------------

options.linAlg = 'standard';
options.linAlg = 'adaptive';
if strcmp(options.linAlg,'standard')
    options.timeStep = 0.05;
    options.taylorTerms = 5;
    options.zonotopeOrder = 50;
elseif strcmp(options.linAlg,'adaptive')
    options.error = 0.05;
end
options.guardIntersect = 'polytope';
options.enclose = {'box'};

reachSet = reach(hybridSys, params, options);


% Simulation --------------------------------------------------------------

simOpt.points = 10;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 10;

simRes = simulateRandom(hybridSys, params, simOpt);


% Visualization -----------------------------------------------------------

% plot different projections
dims = {[1, 2],[1]};

for i = 1:length(dims)

    figure; hold on; box on
    projDims = dims{i};

    % plot reachable sets
    if length(projDims) == 1
        plotOverTime(reachSet, projDims,  'Facecolor', [0.5,0.5,0.5], 'EdgeColor', 'none');
    else
        plot(reachSet, projDims,  'Facecolor', [0.5,0.5,0.5], 'EdgeColor', 'none', 'Filled', true);
    end

    % plot initial set
    if length(projDims) ~= 1
        plot(params.R0, projDims,  'Facecolor', 'white', 'EdgeColor', 'black', 'Filled', true);
    end

    % plot simulation results
    if length(projDims) == 1
        plotOverTime(simRes, projDims, '-', 'color', 'black');
    else
        plot(simRes,projDims, '-', 'color', 'black');
    end

    % label plot
    if length(projDims) == 1
        xlabel('time')
        ylabel(['x_{',num2str(projDims(1)),'}'])
    else
        xlabel(['x_{',num2str(projDims(1)),'}'])
        ylabel(['x_{',num2str(projDims(2)),'}'])
    end

end

end


%% Auxiliary Functions ----------------------------------------------------

function HA = bball(~)


%% Generated on 16-Apr-2021

%---------------Automaton created from Component 'system'------------------

%% Interface Specification:
%   This section clarifies the meaning of state, input & output dimensions
%   by showing their mapping to SpaceEx variable names. 

% Component 1 (system.ball):
%  state x := [x; v]
%  input u := [uDummy]

%-------------------------Component system.ball----------------------------

%-----------------------------State always---------------------------------

%% equation:
%   x' == v & v' == -g
dynA = ...
[0,1;0,0];
dynB = ...
[0;0];
dync = ...
[0;-9.8100000000000004973799150320701];
dynamics = linearSys(dynA, dynB, dync);

%% equation:
%   x >= 0
A = ...
[-1,0];
b = ...
[0];
polyOpt = struct('A', A, 'b', b);
inv = mptPolytope(polyOpt);

trans = {};
%% equation:
%   v' := -c*v
resetA = ...
[1,0;0,-0.75];
resetb = ...
[0;0];
reset = struct('A', resetA, 'b', resetb);

%% equation:
%   x <= eps & v < 0
c = [-1;0];
d = 0;C = ...
[0,1];
D = [0];

guard = conHyperplane(c,d,C,D);

trans{1} = transition(guard, reset, 1);

loc{1} = location('S1', inv, trans, dynamics);



HA = hybridAutomaton(loc);


end
