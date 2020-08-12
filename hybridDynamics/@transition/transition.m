classdef transition
% transition - constructor of the class transition
%
% Syntax:  
%    obj = transition(guard,reset,target)
%
% Inputs:
%    guard - guard set specified as contSet object
%    reset - reset function (only linear map!): Ax+b, with struct:
%            reset.A, reset.b
%    target - target: int (id of target location)
%
% Outputs:
%    obj - generated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton, location

% Author:       Matthias Althoff
% Written:      02-May-2007 
% Last update:  30-July-2016
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    guard = [];         % guard set
    reset = [];         % reset function 
    target = [];        % target location
end

methods
    
    % class constructor
    function obj = transition(varargin)

        % parse input arguments
        if nargin == 3
            obj.guard = varargin{1};
            obj.reset = varargin{2};
            obj.target = varargin{3};
        else
            error('Wrong number of input arguments!');
        end
    end
end
end

%------------- END OF CODE --------------