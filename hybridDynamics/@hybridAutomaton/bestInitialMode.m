function [mode,x0] = bestInitialMode(HA,x,t,varargin)
% bestInitialMode - determines the initial mode that is best suited for 
%   simulating the given trajectory (required because there might exist 
%   multiple suitable modes for one initial state)
%
% Syntax:
%    [mode,x0] = bestInitialMode(HA,x,t)
%    [mode,x0] = bestInitialMode(HA,x,t,u)
%
% Inputs:
%    HA - hybridAutomaton object
%    x - matrix of dimension [n,N] storing the states for the trajectory
%    t - vector of dimension [1,N] storing the time for the trajectory
%    u - matrix of dimension [m,N-1] storing the inputs for the trajectory
%
% Outputs:
%    mode - best initial mode for simulating the automaton
%    x0 - best initial state for simulating the automaton
%
% Example: 
%    HA = roomHeating();
%
%    simOpts.x0 = [20.5; 5];
%    simOpts.tFinal = 2;
%    simOpts.startLoc = 1;
%    [t,x] = simulate(HA,simOpts);
%
%    mode = bestInitialMode(HA,x,t)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: linearSys/identify

% Authors:       Niklas Kochdumper
% Written:       05-June-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    % parse input arguments
    narginchk(3,4);

    if nargin > 3
        u = varargin{1};
    else
        u = [];
    end
    traj = trajectory(u,x,[],t);

    % construct initial state for an ARX model
    x0 = aux_initStateARX(HA,traj);

    % simulation options
    simOpts = struct('x0',x0,'tFinal',traj.t(end));

    if isfield(traj,'u')
        simOpts.u = traj.u;
    end

    % loop over all feasible initial modes
    cost = inf*ones(length(HA.location),1);

    for j = 1:length(HA.location)
        if contains(HA.location(j).invariant,x0)
            cost(j) = aux_simulation(HA.location(j),simOpts,traj);
        end
    end

    % select best simulation
    [~,mode] = min(cost);
end


% Auxiliary functions -----------------------------------------------------

function cost = aux_simulation(sys,simOpts,traj)
% simulate the hybrid automaton and bring data to correct format

    % simulate the hybrid automaton
    [t,x] = simulate(sys.contDynamics,simOpts);

    % bring data to the correct format
    [t,ind] = unique(t);
    x = x(:,ind);

    xInt = interp1(t,x',traj.t,'linear','extrap')';

    % find point where trajectory first leaves the invariant
    ind = size(xInt,2);

    for i = 1:size(xInt,2)
        if ~contains(sys.invariant,xInt(:,i))
            ind = i; break;
        end
    end

    % calculate cost
    n = size(traj.x,1);
    cost = sum(sum((xInt(1:n,1:ind) - traj.x(:,1:ind)).^2,2))/ind;
end

function x0 = aux_initStateARX(HA,traj)
% constrcut the initial state for an ARX model

    % determine number of derivatives for the ARX model
    nAut = dim(HA.location(1).invariant);
    nData = size(traj.x,1);

    ARX = round(nAut/nData) - 1;

    % compute the initial state
    x = traj.x;
    for i = 1:ARX
        dx = aux_derivative(traj.t,x(end-nData+1:end,:));
        x = [x; dx];
    end

    x0 = x(:,1);
end

function [dy,x] = aux_derivative(x,y)
    
    % number of subintervals
    N = length(x)-1;
    
    % preallocates vector to store derivative
    dy = zeros(size(y));
    
    % approximates derivative at lower bound using forward difference
    dy(:,1) = (y(:,2)-y(:,1))/(x(2)-x(1));
    
    % approximates derivative at upper bound using backward difference
    dy(:,N+1) = (y(:,N+1)-y(:,N))/(x(N+1)-x(N));
    
    % approximates derivatives at all other nodes using central differences
    for i = 2:N
        dy(:,i) = (y(:,i+1)-y(:,i-1))/(x(i+1)-x(i-1));
    end  
end

% ------------------------------ END OF CODE ------------------------------
