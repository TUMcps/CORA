classdef location
% location - constructor of class location
%
% Syntax:
%    loc = location()
%    loc = location(invSet,trans,sys)
%    loc = location(name,invSet,trans,sys)
%
% Inputs:
%    name - name of the location
%    invSet - invariant set
%    trans - cell-array containing all transitions
%    sys - continuous dynamics
%
% Outputs:
%    loc - generated location object
%
% Example:
%    % name of location
%    name = 'S1';
%
%    % invariant
%    polyOpt = struct('A',[-1,0],'b',0);
%    inv = mptPolytope(polyOpt);
%    
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%
%    % transition
%    trans{1} = transition(guard,reset,2);
%
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
%
%    % define location
%    loc = location(name,inv,trans,dynamics);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton, transition

% Author:       Matthias Althoff
% Written:      02-May-2007 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    % name of the location
    name (1,:) = 'location';

    % invariant set 
    invariant = [];

    % cell-array storing the transitions
    transition (:,1) cell = {};

    % system dynamics
    contDynamics {mustBeA(contDynamics,'contDynamics')} = contDynamics();
end

methods
    
    % class constructor
    function loc = location(varargin)

        % parse input arguments
        if nargin == 0
            invSet = [];
        elseif nargin <= 2
            throw(CORAerror('CORA:notEnoughInputArgs',3));
        elseif nargin == 3
            invSet = varargin{1};
            loc.transition = varargin{2};
            loc.contDynamics = varargin{3};
        elseif nargin == 4
            loc.name = varargin{1};
            invSet = varargin{2};
            loc.transition = varargin{3};
            loc.contDynamics = varargin{4}; 
        else
            throw(CORAerror('CORA:tooManyInputArgs',4));
        end

        % elements of transition cell-array have to be transition objects
        if ~all(cellfun(@(x) isa(x,'transition'),loc.transition,'UniformOutput',true))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'All elements of transition have to be transition objects.'));
        end

        % convert invariant sets to polytopes if possible
        if (isnumeric(invSet) && isempty(invSet)) || ...
                isa(invSet,'mptPolytope') || isa(invSet,'levelSet')
            % keep mptPolytopes, levelSet, and empty invariants as is
            loc.invariant = invSet;
        else
            loc.invariant = mptPolytope(invSet);
        end
    end
end
end

%------------- END OF CODE --------------