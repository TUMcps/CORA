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
%    trans - object-array containing all transitions
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
%    inv = polytope([-1,0],0);
%    
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%
%    % reset function
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%
%    % transition
%    trans(1) = transition(guard,reset,2);
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

% Authors:       Matthias Althoff
% Written:       02-May-2007 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    % name of the location
    name (1,:) = 'location';

    % invariant set 
    invariant = [];

    % array storing the transitions
    transition (:,1) = transition();

    % system dynamics
    contDynamics = contDynamics();
end

methods
    
    % class constructor
    function loc = location(varargin)

        % parse input arguments
        if nargin == 0
            return
        elseif nargin == 1 && isa(varargin{1},'location')
            % copy constructor
            loc = varargin{1}; return
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

        % loc.transition has to be an array of transition objects
        if ~isa(loc.transition,'transition')
            if iscell(loc.transition)
                % legacy error message
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Transitions have to be an object array instead of a cell array.'));
            else
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Transitions have to be a transition object array.'));
            end
        end

        % convert invariant sets to polytopes if possible
        if isa(invSet,'fullspace') || isa(invSet,'polytope') ...
                || isa(invSet,'levelSet') || isa(invSet,'emptySet')
            % keep polytopes, levelSet, and empty invariants as is
            loc.invariant = invSet;
        else
            loc.invariant = polytope(invSet);
        end
    end
end
end

% ------------------------------ END OF CODE ------------------------------
