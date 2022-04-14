classdef location
% location - constructor of class location
%
% Syntax:  
%    obj = location(invSet,trans,sys)
%    obj = location(name,invSet,trans,sys)
%
% Inputs:
%    name - name of the location: char array
%    invSet - invariant set (class: contSet)
%    trans - cell-array containing all transitions
%    sys - continuous dynamics (class: contDynamics) 
%
% Outputs:
%    obj - gnerated object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton, transition

% Author: Matthias Althoff
% Written: 02-May-2007 
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    name = [];              % name of the location
    invariant = [];         % invariant set 
    transition = [];        % cell-array storing the transitions
    contDynamics = [];      % system dynamics (contDynamics object)
end

methods
    
    % class constructor
    function obj = location(varargin)

        % parse input arguments
        if nargin == 3
           name = 'location';
           invSet = varargin{1};
           trans = varargin{2};
           sys = varargin{3};
        elseif nargin == 4
           name = varargin{1};
           invSet = varargin{2};
           trans = varargin{3};
           sys = varargin{4}; 
        else
            error('Wrong number of inputs arguments!');
        end

        % assign object properties
        obj.name = name;
        obj.transition = trans;
        obj.contDynamics = sys;

        % convert invariant sets to polytopes if possible
        if isa(invSet,'mptPolytope') || isa(invSet,'levelSet')
            obj.invariant = invSet;
        else
            obj.invariant = mptPolytope(invSet);
        end
    end
end
end

%------------- END OF CODE --------------