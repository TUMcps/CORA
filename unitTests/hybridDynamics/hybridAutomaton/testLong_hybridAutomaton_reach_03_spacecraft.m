function res = testLong_hybridAutomaton_reach_03_spacecraft
% testLong_hybridAutomaton_reach_03_spacecraft - unit test function
%    for hybrid systems  with nonlinear dynamics based on the spacecraft
%    rendezvous example [1]
%       
% Syntax:
%    res = testLong_hybridAutomaton_reach_03_spacecraft
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% References: 
%   [1] N. Chan et al. "Verifying safety of an autonomous spacecraft 
%       rendezvous mission (Benchmark proposal)"  

% Authors:       Niklas Kochdumper
% Written:       23-December-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Reachability Settings ---------------------------------------------------

% problem description
R0 = zonotope([[-900;-400;0;0;0],diag([25;25;0;0;0])]);

params.R0 = R0;                                     % initial set
params.startLoc = 1;                                % initial location
params.finalLoc = 0;                                % 0: no final location
params.tFinal = 300;                                % final time

% uncertain input
params.U = zonotope(0); 


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.timeStep = 2e-1;
options.taylorTerms = 30;
options.zonotopeOrder = 10;
options.reductionTechnique = 'girard';

% copy some options (since we call location/reach)
options.U = params.U;
options.R0 = params.R0;
options.finalLoc = params.finalLoc;
options.tFinal = params.tFinal;
options.tStart = 0;
options.startLoc = params.startLoc;
options.specification = [];
options.intersectInvariant = false;


% Hybrid Automaton --------------------------------------------------------

% specify hybrid automata
HA = rendezvous_SRNA01();           % automatically converted from SpaceEx

% get first location
loc = HA.location(1);


% Simulation --------------------------------------------------------------

% settings for random simulation
simOpt.points = 4;         % number of initial points
simOpt.fracVert = 0.5;     % fraction of vertices initial set
simOpt.fracInpVert = 0.5;  % fraction of vertices input set
simOpt.nrConstInp = 2;     % number of constant input over time horizon  

% random simulation
simRes = simulateRandom(HA,params,simOpt); 

% extract the states at the intersection with the guard set
[~,xHit] = extractHits(simRes,1);
points = reshape(cell2mat(xHit),5,[]);
% project for verification
points = points(2:end,:);


% Reachability Analysis ---------------------------------------------------

guardIntersect = {'polytope','zonoGirard','conZonotope', ...
                  'hyperplaneMap','pancake','nondetGuard'};

res = true(length(guardIntersect),1);
for i = 1:length(guardIntersect)

    % set options for guard intersection 
    options_ = options;
    options_.guardIntersect = guardIntersect{i};
    
    switch guardIntersect{i}
        
        case 'polytope'
            options_.enclose = {'pca'};
            
        case 'zonoGirard'
            options_.enclose = {'pca'};
            
        case 'conZonotope'
            options_.enclose = {'pca'};
            options_.guardOrder = 5;
            
        case 'hyperplaneMap'
            options_.guardOrder = 3;
            
        case 'nondetGuard'
            options_.enclose = {'pca'};
    end
    
    % reachability analysis    
    [~,Rjump,~] = reach(loc,params.R0,interval(options.tStart),options_);
    
    % check if all simulated points are located inside the computed guard
    % intersection
    P = polytope(project(Rjump(1).set,2:5));
    res(i) = all(contains(P,points,'exact',1e-3));
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
