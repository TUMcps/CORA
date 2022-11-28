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
% Last revision:---

%------------- BEGIN CODE --------------
  

properties (SetAccess = private, GetAccess = public)
    % name of the system
    name (1,:) = '';

    % state dimension
    dim (1,1) {mustBeInteger,mustBeNonnegative} = 0;

    % input dimension
    nrOfInputs (1,1) {mustBeInteger,mustBeNonnegative} = 0;

    % output dimension
    nrOfOutputs (1,1) {mustBeInteger,mustBeNonnegative} = 0;
end
    
methods
    
    % class constructor
    function obj = contDynamics(varargin)
        
        if nargin == 0
            % empty
        elseif nargin == 1
            obj.name = varargin{1};
        elseif nargin == 2
            obj.name = varargin{1};
            obj.dim = varargin{2};
        elseif nargin == 3
            obj.name = varargin{1};
            obj.dim = varargin{2};
            obj.nrOfInputs = varargin{3};
        elseif nargin == 4
            obj.name = varargin{1};
            obj.dim = varargin{2};
            obj.nrOfInputs = varargin{3};
            obj.nrOfOutputs = varargin{4};
        else
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end
    end
end

end

%------------- END OF CODE --------------