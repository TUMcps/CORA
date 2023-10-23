function res = test_spaceex2cora_hybrid_parallel_01_emptycomp
% test_spaceex2cora_hybrid_parallel_01_emptycomp - test for model
%    conversion from SpaceEx to CORA for a parallel hybrid automaton with
%    synchronization labels, where one component does not have any flow
%    equation or invariant set, but transitions with synchronization labels
%
% Syntax:
%    test_spaceex2cora_hybrid_parallel_01_emptycomp
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       15-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;


% directory to SpaceEx model file
dir_spaceex = [CORAROOT filesep 'unitTests' filesep 'converter' ...
    filesep 'spaceex2cora' filesep 'testSystems'];

% file name of SpaceEx model file
filename = 'test_hybrid_parallel_emptycomp';

% convert SpaceEx model from .xml file
spaceex2cora([dir_spaceex filesep filename '.xml']);

% instantiate system from converted SpaceEx model
sys_spaceex = feval(filename);


% instantiate equivalent CORA model

% first component, first location
dynamics = linearSys([0 0; 0 0],[0; 0],[0; 1]);
inv = polytope([0 1],0);
guard = conHyperplane([0 1],0);
reset = struct('A',[1 0; 0 1],'c',[1;5]);
trans = transition(guard,reset,2,'switch12');
loc(1) = location('loc1',inv,trans,dynamics);

% first component, second location
dynamics = linearSys([0 0; 0 0],[0; 0],[0; -1]);
inv = polytope([0 -1],0);
guard = conHyperplane([0 1],0);
reset = struct('A',[1 0; 0 1],'c',[1;-5]);
trans = transition(guard,reset,1,'switch21');
loc(2) = location('loc2',inv,trans,dynamics);

% init first hybrid automaton
HA1 = hybridAutomaton(loc);

% second component, first location
dynamics = linearSys([],zeros(0,1));
inv = fullspace(0);
guard = fullspace(0);
reset = struct('A',[],'c',zeros(0,1));
trans = transition(guard,reset,2,'switch12');
loc(1) = location('loc1',inv,trans,dynamics);

% second component, second location
dynamics = linearSys([],zeros(0,1));
inv = fullspace(0);
guard = fullspace(0);
reset = struct('A',[],'c',zeros(0,1));
trans = transition(guard,reset,1,'switch21');
loc(2) = location('loc2',inv,trans,dynamics);

% init second hybrid automaton
HA2 = hybridAutomaton(loc);

% components and input binds
components = [HA1;HA2];
inputBinds{1} = [0 1];
inputBinds{2} = [0 1];

% instantiate parallel hybrid automaton
sys_cora = parallelHybridAutomaton(components,inputBinds);

% compare systems
if ~isequal(sys_cora,sys_spaceex,1e-4)
    res = false;
end


% execution of the model

% model parameters
params.R0 = zonotope([-2;-3],[0.1; 0]);
params.U = zonotope(0);
params.startLoc = [1;1];
params.tFinal = 20;

% settings for continuous reachability 
options.timeStep = 1;
options.taylorTerms = 4;
options.zonotopeOrder = 5;
options.verbose = true;

% settings for hybrid systems
options.guardIntersect = 'polytope';
options.enclose = {'box'};

% reachability analysis
R = reach(sys_spaceex,params,options);

% final set
R_end = query(R,'finalSet');
R_end = R_end{1};

% final set according to manual computation
R_end_ = zonotope([2;-3],[0.1; 0]);

if ~isequal(R_end,R_end_)
    res = false;
end

% visualization (debugging)
figure; hold on; box on;
plot(R,[1,2]);
plot(params.R0,[1,2],'FaceColor',colorblind('gray'),'EdgeColor','k');
plot(R_end_,[1,2],'g','LineWidth',1.5);
plot(R_end,[1,2],'r--','LineWidth',1.5);
close;

% ------------------------------ END OF CODE ------------------------------
