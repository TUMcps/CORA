function res = simulateConstrainedRandom(obj, options)
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
%    res = simulateConstrainedRandom(obj, options)
%
% Inputs:
%    obj - contDynamics object
%    options - model parameters and settings for random simulation
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.

% Authors:       Matthias Althoff
% Written:       05-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% initialize time and state
t = cell(options.points,1);
x = cell(options.points,1);

% output equation only for linearSys and linearSysDT currently
comp_y = (isa(obj,'linearSys') || isa(obj,'linearSysDT')) && ~isempty(obj.C);
if comp_y; y = cell(options.points,1); end

% generate random initial points
nrExtreme = ceil(options.points*options.fracVert);
nrStandard = options.points - nrExtreme;
X0 = [];
if nrExtreme > 0
	X0 = [X0, randPoint(options.R0,nrExtreme,'extreme')]; 
end
if nrStandard > 0
	X0 = [X0, randPoint(options.R0,nrStandard,'standard')];
end

% loop over all starting points in X0
for r = 1:options.points
    
    % Is output desired?
    if comp_y
        y{r} = zeros(1,obj.nrOfOutputs); 
    end
    
    % start of trajectory
    options.x0 = X0(:,r);

    % simulate dynamical system
    if comp_y
        [t{r},x{r},~,y{r}] = simulateConstrained(obj,options);
    else
        [t{r},x{r}] = simulateConstrained(obj,options);
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
