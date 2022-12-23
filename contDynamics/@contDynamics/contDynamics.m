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

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      02-May-2007 
% Last update:  18-March-2016
%               04-August-2016 (changed to new OO format)
%               04-March-2019 (number of outputs added)
%               22-May-2020 (NK, deleted stateID, inputID, etc. properties)
%               14-December-2022 (TL, property check in inputArgsCheck)
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    % name of the system
    name;

    % state dimension
    dim;

    % input dimension
    nrOfInputs;

    % output dimension
    nrOfOutputs;
end
    
methods
    
    % class constructor
    function obj = contDynamics(varargin)

        % parse input
        maxArgs = 4; % max input args
        if nargin > maxArgs
            throw(CORAerror('CORA:tooManyInputArgs',maxArgs));
        end
        [name, dim, nrOfInputs, nrOfOutputs] = setDefaultValues( ...
            {'',0,0,0}, varargin);
        inputArgsCheck({ ...
            {name, 'att', {'char', 'string'}}; ...
            {dim, 'att', 'numeric', ...
                {'integer', 'nonnegative', 'scalar'}}; ...
            {nrOfInputs, 'att', 'numeric', ...
                {'integer', 'nonnegative', 'scalar'}}; ...
            {nrOfOutputs, 'att', 'numeric', ...
                {'integer', 'nonnegative', 'scalar'}}; ...
        })

        % set properties
        obj.name = name;
        obj.dim = dim;
        obj.nrOfInputs = nrOfInputs;
        obj.nrOfOutputs = nrOfOutputs;
    end
end

end

%------------- END OF CODE --------------