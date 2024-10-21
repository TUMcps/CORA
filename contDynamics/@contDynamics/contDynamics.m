classdef contDynamics < matlab.mixin.Copyable
% contDynamics - basic class for continuous dynamics
%
% Syntax:
%    sys = contDynamics()
%    sys = contDynamics(name)
%    sys = contDynamics(name,states)
%    sys = contDynamics(name,states,inputs)
%    sys = contDynamics(name,states,inputs,outputs)
%    sys = contDynamics(name,states,inputs,outputs,dists)
%    sys = contDynamics(name,states,inputs,outputs,dists,noises)
%
% Inputs:
%    name - system name
%    states - number of states
%    inputs - number of inputs
%    outputs - number of outputs
%    dists - number of disturbances
%    noises - number of disturbances on output
%
% Outputs:
%    sys - generated contDynamics object
%
% Example:
%    sys = contDynamics('system',2,1,1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       02-May-2007 
% Last update:   18-March-2016
%                04-August-2016 (changed to new OO format)
%                04-March-2019 (number of outputs added)
%                22-May-2020 (NK, deleted stateID, inputID, etc. properties)
%                14-December-2022 (TL, property check in inputArgsCheck)
%                30-August-2024 (MW, add disturbance/noise sizes)
%                30-August-2024 (MW, add disturbance/noise sizes)
%                16-October-2024 (TL, renames dim to nrOfStates)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    name;               % name of the system
    nrOfStates;         % state dimension
    nrOfInputs;         % input dimension
    nrOfOutputs;        % output dimension
    nrOfDisturbances;   % disturbance dimension
    nrOfNoises;         % noise dimension

    % legacy
    dim;                % (also) state dimension
end
    
methods
    
    % class constructor
    function sys = contDynamics(varargin)

        % 0. check number of input arguments
        assertNarginConstructor(0:6,nargin);

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'contDynamics')
            sys = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [name,states,inputs,outputs,dists,noises] = ...
            setDefaultValues({'',0,0,0,0,0},varargin);

        % 3. check correctness of input arguments
        if CHECKS_ENABLED && nargin > 0
            inputArgsCheck({ ...
                {name, 'att', {'char', 'string'}}; ...
                {states, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {inputs, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {outputs, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {dists, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {noises, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
            })
        end

        % 4. assign properties
        sys.name = name;
        sys.nrOfStates = states;
        sys.nrOfInputs = inputs;
        sys.nrOfOutputs = outputs;
        sys.nrOfDisturbances = dists;
        sys.nrOfNoises = noises;
    end
end

methods (Access = protected)
    [printOrder] = getPrintSystemInfo(S)
end

% getter & setter ---------------------------------------------------------

methods 
    function dim = get.dim(sys)
        CORAwarning('CORA:deprecated', 'property', 'contDynamics.dim', 'CORA v2025', ...
            'Please use contDynamics.nrOfStates instead.', ...
            'This change was made to be consistent with the other properties.')
        dim = sys.nrOfStates;
    end
    function set.dim(sys,dim)
        CORAwarning('CORA:deprecated', 'property', 'contDynamics.dim', 'CORA v2025', ...
            'Please use contDynamics.nrOfStates instead.', ...
            'This change was made to be consistent with the other properties.')
        sys.nrOfStates = dim;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
