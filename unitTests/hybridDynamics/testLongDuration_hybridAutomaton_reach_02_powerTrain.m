function res = testLongDuration_hybridAutomaton_reach_02_powerTrain
% testLongDuration_hybridAutomaton_reach_02_powerTrain - unit test function for
%    hybrid systems with linear continuous dynamics
%
% Checks the solution of the hybrid system class for the powertrain example
% in [1]. It is checked if the computed guard intersection encloses all
% simulated trajectories.
%
% Syntax:  
%    res = testLongDuration_hybridAutomaton_reach_02_powerTrain
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% References:
%   [1] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems"

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-September-2011
% Last update:  23-December-2019
% Last revision:---

%------------- BEGIN CODE --------------


% Hybrid Automaton --------------------------------------------------------

% load power train parameters
[HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(11); 

% get the relevant location
loc = HA.location{3};


% Parameter ---------------------------------------------------------------

% problem description
params.R0 = zonotope([Zcenter, Zdelta]);       % initial set
params.startLoc = 3;                           % initial location
params.finalLoc = 0;                           % 0: no final location


% Reachability Settings ---------------------------------------------------

% settings for continuous reachability 
options.timeStep = 5e-4; 
options.taylorTerms = 20;
options.zonotopeOrder = 20;


% Simulation (First Part) -------------------------------------------------

% system input:
u = [-5; 0];                                % acceleration and load torque

simOptions.u{1} = B1*u + c1; 
simOptions.u{2} = B2*u + c2; 
simOptions.u{3} = B3*u + c3; 

% set start and end time
params.tFinal = 0.2;
simOptions.tFinal = params.tFinal;
simOptions.startLoc = 3;

% obtain random simulation results
N = 10;
simPoints = cell(N,1);

for i = 1:N
    
    % set initial state
    if i == 1
        simOptions.x0 = Zcenter + Zdelta;
    elseif i == 2
        simOptions.x0 = Zcenter - Zdelta;
    else
        simOptions.x0 = randPoint(params.R0);
    end 
    
    % simulate hybrid automaton
    [t,x,locID] = simulate(HA,simOptions); 
    
    % store results
    simPoints{i} = x{end}(end,:)';
    
    if i == 1
        simRes1 = simResult(x,t,locID); 
    else
        simRes1 = add(simRes1,simResult(x,t,locID)); 
    end
end


% Simulation (Second Part) ------------------------------------------------

% system input:
u = [5; 0];                                 % acceleration and load torque

simOptions.u{1} = B1*u + c1; 
simOptions.u{2} = B2*u + c2; 
simOptions.u{3} = B3*u + c3; 

% set start and final time
params.tStart = 0.2;
params.tFinal = 2;
simOptions.tStart = params.tStart;
simOptions.tFinal = params.tFinal;

% obtain random simulation results
for i = 1:N
    
    % set initial state
    simOptions.x0 = simPoints{i};
    
    % simulate hybrid automaton
    [t,x,locID] = simulate(HA,simOptions); 
    
    % store resulets
    if i == 1
        simRes2 = simResult(x,t,locID); 
    else
        simRes2 = add(simRes2,simResult(x,t,locID)); 
    end
end

% extract the states at the intersection with the guard set
points = [];

for j = 1:length(simRes2.x)
   if simRes2.loc{j} == 3
       p = simRes2.x{j}(end,:)';
       points = [points,p];
   end
end


% Reachability Analysis ---------------------------------------------------

guardIntersect = {'polytope','zonoGirard','conZonotope', ...
                  'hyperplaneMap','pancake','nondetGuard'};

for i = 1:length(guardIntersect)

    % set options for guard intersection 
    options_ = options;
    params_ = params;
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
    
    % reachability analysis (first part)
    u = [-5; 0];                                

    params_.U{1} = zonotope(B1*u + c1); 
    params_.U{2} = zonotope(B2*u + c2); 
    params_.U{3} = zonotope(B3*u + c3); 

    params_.tStart = 0;
    params_.tFinal = 0.2;
    
    R = reach(HA,params_,options_);

    % reachability analysis (second part)
    u = [5; 0];                                

    options_.U = zonotope(B3*u + c3);
    options_.R0 = R.timePoint.set{end};
    options_.finalLoc = params.finalLoc;
    options_.startLoc = params.startLoc;
    options_.tFinal = 2;
    options_.tStart = 0.2;
    options_.reductionTechnique = 'girard';
    options_.specification = [];
    
    [~,Rjump,~] = reach(loc,options_.R0,interval(options_.tStart),options_);
    
    % check if all simulated points are located inside the computed guard
    % intersection
    R = project(Rjump{1}.set,2:11);
    if size(R.Z,2) > 12
        R = reduce(R,'pca',1);
    end
    poly = mptPolytope(R);
    
    points_ = points(2:end,:);
    
    for j = 1:size(points_,2)
       if ~in(poly,points_(:,j),1e-10)
          res = 0;
          error('Guard intersection with method %s failed!',guardIntersect{i});
       end
    end
end

res = 1;

%------------- END OF CODE --------------