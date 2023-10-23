function res = example_hybrid_reach_02_powerTrain()
% example_hybrid_reach_02_powerTrain - power train example from [1]
%
% Syntax:
%    res = example_hybrid_reach_02_powerTrain()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% References:
%   [1] M. Althoff et al. "Avoiding Geometic Intersection Operations in 
%       Reachability Analysis of Hybrid Systems"

% Authors:       Matthias Althoff
% Written:       21-September-2011
% Last update:   23-December-2019
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------


% Hybrid Automaton --------------------------------------------------------

[HA, Zcenter, Zdelta, B1, B2, B3, c1, c2, c3] = initPowerTrain(11); 


% Parameter ---------------------------------------------------------------

% problem description
params.R0 = zonotope([Zcenter, Zdelta]);       % initial set
params.startLoc = 3;                           % initial location


% Reachability Options ----------------------------------------------------

% settings for continuous reachability 
options.timeStep = 5e-4; 
options.taylorTerms = 20;
options.zonotopeOrder = 20;

% settings for hybrid systems
options.guardIntersect = 'zonoGirard';
options.enclose = {'pca'};


% Simulation Options ------------------------------------------------------

simOptions.startLoc = params.startLoc; 


% FIRST PART --------------------------------------------------------------

% system input:
u = [-5; 0];                                % acceleration and load torque

params.U{1} = zonotope(B1*u + c1);
simOptions.u{1} = B1*u + c1; 

params.U{2} = zonotope(B2*u + c2);
simOptions.u{2} = B2*u + c2; 

params.U{3} = zonotope(B3*u + c3);
simOptions.u{3} = B3*u + c3; 

% set start and end time
params.tFinal = 0.2;
simOptions.tFinal = params.tFinal;

% obtain random simulation results
N = 10;
simPoints = cell(N,1);
simRes = [];

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
    [t,x,loc] = simulate(HA,simOptions); 
    
    % store results
    simPoints{i} = x{end}(end,:)';
    simRes = add(simRes,simResult(x,t,loc)); 
end

% compute reachable set
tStart = tic;
R = reach(HA,params,options);
tComp = toc(tStart);

disp(['Computation time (first part): ',num2str(tComp),'s']);


% SECOND PART -------------------------------------------------------------

% system input:
u = [5; 0];                                 % acceleration and load torque

params.U{1} = zonotope(B1*u + c1);
simOptions.u{1} = B1*u + c1; 

params.U{2} = zonotope(B2*u + c2);
simOptions.u{2} = B2*u + c2; 

params.U{3} = zonotope(B3*u + c3);
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
    [t,x,loc] = simulate(HA,simOptions); 
    
    % store resulets
    simRes = add(simRes,simResult(x,t,loc));
end

% compute reachable set
params.R0 = R(end).timePoint.set{end};

tStart = tic;
Rtemp = reach(HA,params,options);
R = add(R,Rtemp);
tComp = toc(tStart);

disp(['Computation time (second part): ',num2str(tComp),'s']);


% Visualization -----------------------------------------------------------

dims = {[1,2],[1,3]};

for i = 1:length(dims)
    
    figure; hold on; box on;
    useCORAcolors("CORA:contDynamics")
    projDim = dims{i};
    
    % plot reachable set
    plot(R,projDim,'Order',3);
    
    updateColorIndex

    % plot simulation results
    plot(simRes,projDim);

    plot([-0.03 -0.03],[-40,100],'--r','DisplayName','Guard set')
    plot([0.03 0.03],[-40,100],'--r','HandleVisibility','off')
    
    % formatting
    xlabel(['x_',num2str(projDim(1))]);
    ylabel(['x_',num2str(projDim(2))]);
end

% example completed
res = true;

% ------------------------------ END OF CODE ------------------------------
