function [res, varargout] = simulateRandom(obj, params, options,varargin)
% simulateRandom - performs several random simulation of the system. It 
% can be set how many simulations should be performed, what percentage of 
% initial states should start at vertices of the initial set, what 
% percentage of inputs should be chosen from vertices of the input set, and 
% how often the input should be changed.
%
% Syntax:  
%   res = simulateRandom(obj, params, options)
%   res = simulateRandom(obj, params, options, type)
%   res = simulateRandom(obj, params, options, 'gaussian', R)
%
% Inputs:
%    obj - contDynamics object
%    params - system parameters
%    options - settings for random simulation, depending on type (see below)
%       type = 'gaussian' or 'standard':
%           .points - nr of simulation runs
%           .fracVert - fraction of initial states starting from vertices
%           .fracInpVert - fraction of input values taken from the 
%                          vertices of the input set
%           .inpChanges - number of times the input is changed in a simulation run
%       type = 'rrt':
%           .points:    number of random initial points (positive integer)
%           .vertSamp:  flag that specifies if random initial points, inputs,
%                       and parameters are sampled from the vertices of the 
%                       corresponding sets (0 or 1)
%           .strechFac: stretching factor for enlarging the reachable sets 
%                       during execution of the algorithm (scalar > 1).
%    type - type of random sample generation ('standard', 'rrt' and 'gaussian')
%
%    for type = 'rrt', additional input argument:
%    R - object of class reachSet storing the computed reachable set
%
% Outputs:
%    res - object of class simResult storing time and states of the 
%          simulated trajectories.
%
%    for type = 'rrt', additional output argument:
%    points - final points of the simulation
%
% 
% Author:       Matthias Althoff
% Written:      17-August-2016
% Last update:  08-May-2020 (MW, update interface)
%               28-June-2021 (MP, unify random simulation functions)
% Last revision:---


%------------- BEGIN CODE --------------

% set default type option
type = 'standard';

% parse inputs
if nargin > 3 && ~isempty(varargin{1})
    if isa(varargin{1},'char')
        type = varargin{1};
    else
        error(errWrongInput('type'));
    end
end

if nargin > 4 && ~isempty(varargin{2}) 
    if strcmp(type,'rrt')
        R = varargin{2};
    else
        error('Types ''standard'' and ''gaussian'' only support up to 4 input arguments!');
    end
end

% call private simulation function based on type
if strcmp(type,'standard')
    res = simulateStandard(obj,params,options);
elseif strcmp (type,'gaussian')
    res = simulateGaussian(obj,params,options);
elseif strcmp (type,'rrt')
    [res,varargout{1}] = simulateRRT(obj,R,params,options);
else
    error('Incorrect value for argument ''type''!');
end
%------------- END OF CODE --------------