function res = simulateConstrainedRandom(sys,params,options)
% simulateConstrainedRandom - performs several random simulation of the system
%    so that the simulations stay within a given reachable set; this 
%    reachable set is typically a backwards minmax reachable set. These 
%    reachable sets assume that a control input exists such that the 
%    solutions stay within the reachable set. However, the inputs are yet
%    to be determined by simulateConstrained of the individual system 
%    dynamics classes. It can be set how many simulations should be 
%    performed and what percentage of initial states should start at 
%    vertices of the initial set.
%
% Syntax:
%    res = simulateConstrainedRandom(sys,params,options)
%
% Inputs:
%    sys - contDynamics object
%    params - model parameters
%    options - settings for random simulation
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff
% Written:       05-January-2023
% Last update:   29-June-2024 (TL, bug fix for empty R0)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate random initial points
nrExtreme = ceil(options.points*options.fracVert);
nrStandard = options.points - nrExtreme;
X0 = [];
if nrExtreme > 0
	X0 = [X0, randPoint(params.R0,nrExtreme,'extreme')]; 
end
if nrStandard > 0
	X0 = [X0, randPoint(params.R0,nrStandard,'standard')];
end

% check number of generated points (might be different for e.g. empty sets)
nrPoints = size(X0,2);

% initialize time and state
t = cell(nrPoints,1);
x = cell(nrPoints,1);
% output equation only for linearSys and linearSysDT currently
comp_y = (isa(sys,'linearSys') || isa(sys,'linearSysDT')) && ~isempty(sys.C);
if comp_y; y = cell(nrPoints,1); end

% loop over all starting points in X0
for r = 1:nrPoints
    
    % Is output desired?
    if comp_y
        y{r} = zeros(1,sys.nrOfOutputs); 
    end
    
    % start of trajectory
    params.x0 = X0(:,r);

    % simulate dynamical system
    if comp_y
        [t{r},x{r},~,y{r}] = simulateConstrained(sys,params,options);
    else
        [t{r},x{r}] = simulateConstrained(sys,params,options);
    end
    
end

% construct object storing the simulation results
if comp_y
    res = simResult(x,t,{},y);
else
    res = simResult(x,t);
end

end

% ------------------------------ END OF CODE ------------------------------
