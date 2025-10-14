function traj = simulateRandom(sys, params, varargin)
% simulateRandom - performs several random simulation of the system. It 
%    can be set how many simulations should be performed, what percentage
%    of initial states should start at vertices of the initial set, what 
%    percentage of inputs should be chosen from vertices of the input set,
%    and of how many piecewise constant parts the input is constructed.
%
% Syntax:
%    res = simulateRandom(sys, params, options)
%
% Inputs:
%    sys - contDynamics object
%    params - system parameters
%    options - settings for random simulation, depending on .type (see below)
%       .type = 'gaussian', 'standard' (default), 'rrt', 'constrained';
%       .points - nr of simulation runs
%       further options if .type = 'standard':
%           .fracVert - fraction of initial states starting from vertices
%           .fracInpVert - fraction of input values taken from the 
%                          vertices of the input set
%           .nrConstInp - number of piecewise-constant input segments
%       further options if .type = 'gaussian':
%           .p_conf - probability that a value is within the set
%           .nrConstInp - number of piecewise-constant input segments
%       further options if .type = 'rrt':
%           .points:     number of random initial points (positive integer)
%           .vertSamp:   flag that specifies if random initial points, inputs,
%                        and parameters are sampled from the vertices of the 
%                        corresponding sets (0 or 1)
%           .stretchFac: stretching factor for enlarging the reachable sets 
%                        during execution of the algorithm (scalar > 1).
%           .R:          object of class reachSet storing the computed reachable set
%       further options if .type = 'constrained':
%           .R:          object of class reachSet storing the computed reachable set
%
% Outputs:
%    traj - object of class trajectory storing time and states of the 
%             simulated trajectories
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff
% Written:       17-August-2016
% Last update:   08-May-2020 (MW, update interface)
%                28-June-2021 (MP, unify random simulation functions)
%                05-January-2023 (MA, constrained simulation added)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% input argument validation
options = setDefaultValues({struct()},varargin);
[params,options] = validateOptions(sys,params,options);

% call private simulation function based on type
switch options.type
    case 'standard'
        traj = priv_simulateStandard(sys,params,options);
    case 'gaussian'
        traj = priv_simulateGaussian(sys,params,options);
    case 'rrt'
        traj = priv_simulateRRT(sys,params,options);
    case 'constrained'
        traj = priv_simulateConstrainedRandom(sys,params,options);
end

% ------------------------------ END OF CODE ------------------------------
