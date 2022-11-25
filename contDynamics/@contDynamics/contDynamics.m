classdef contDynamics < matlab.mixin.Copyable
% contDynamics class - basic class for continuous dynamics
%
% Syntax:  
%    obj = contDynamics(name,states,inputs,outputs)
%
% Inputs:
%    name - name of the continuous dynamics: char array
%    states - number of states
%    inputs - number of inputs
%    outputs - number of outputs
%
% Outputs:
%    obj - Generated Object
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
    name = [];
    dim = [];
    nrOfInputs = [];
    nrOfOutputs = [];
end
    
methods
    
    % class constructor
    function obj = contDynamics(varargin)
        
        if nargin==1
            obj.name = varargin{1};
        elseif nargin==4
            obj.name = varargin{1};
            obj.dim = varargin{2};
            obj.nrOfInputs = varargin{3};
            obj.nrOfOutputs = varargin{4};
        else
            error('Wrong number of input arguments');
        end
    end
end

end

%------------- END OF CODE --------------