classdef contDynamics < matlab.mixin.Copyable
% contDynamics - basic class for continuous dynamics
%
% Syntax:
%    obj = contDynamics()
%    obj = contDynamics(name)
%    obj = contDynamics(name,states)
%    obj = contDynamics(name,states,inputs)
%    obj = contDynamics(name,states,inputs,outputs)
%
% Inputs:
%    name - system name
%    states - number of states
%    inputs - number of inputs
%    outputs - number of outputs
%
% Outputs:
%    obj - generated contDynamics object
%
% Example:
%    obj = contDynamics('system',2,1,1);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       02-May-2007 
% Last update:   18-March-2016
%                04-August-2016 (changed to new OO format)
%                04-March-2019 (number of outputs added)
%                22-May-2020 (NK, deleted stateID, inputID, etc. properties)
%                14-December-2022 (TL, property check in inputArgsCheck)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
  

properties (SetAccess = private, GetAccess = public)
    name;           % name of the system
    dim;            % state dimension
    nrOfInputs;     % input dimension
    nrOfOutputs;    % output dimension
end
    
methods
    
    % class constructor
    function obj = contDynamics(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'contDynamics')
            obj = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        if nargin > 4
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end
        [name,n,m,r] = setDefaultValues({'',0,0,0},varargin);

        % 3. check correctness of input arguments
        if CHECKS_ENABLED && nargin > 0
            inputArgsCheck({ ...
                {name, 'att', {'char', 'string'}}; ...
                {n, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {m, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
                {r, 'att', 'numeric', ...
                    {'integer', 'nonnegative', 'scalar'}}; ...
            })
        end

        % 4. assign properties
        obj.name = name;
        obj.dim = n;
        obj.nrOfInputs = m;
        obj.nrOfOutputs = r;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
