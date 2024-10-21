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
%    % transition: guard set, reset function, target
%    guard = polytope([0,1],0,[-1,0],0);
%    reset = linearReset([1,0;0,-0.75]);
%    trans = transition(guard,reset,2);
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
% Last revision: 14-October-2024 (MW, update to current constructor structure)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)

    name = 'location';	            % name of the location
    invariant = [];                 % invariant set 
    transition = transition();      % array storing the transitions
    contDynamics = contDynamics();  % system dynamics

end

methods
    
    % class constructor
    function loc = location(varargin)

        % 0. empty
        assertNarginConstructor([0,1,3,4],nargin);
        if nargin == 0
            return
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'location')
            loc = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [name,inv,trans,sys] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,inv,trans,sys,nargin);

        % 4. compute dependent properties
        [name,inv,trans,sys] = aux_computeProperties(name,inv,trans,sys);

        % 5. assign properties
        loc.name = name;
        loc.invariant = inv;
        loc.transition = trans;
        loc.contDynamics = sys;
        
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [name,inv,trans,sys] = aux_parseInputArgs(varargin)

    % default properties
    name = 'location'; inv = []; trans = transition(); sys = contDynamics();

    % parse arguments
    if nargin == 3
        [inv,trans,sys] = setDefaultValues({inv,trans,sys},varargin);
    elseif nargin == 4
        [name,inv,trans,sys] = setDefaultValues({name,inv,trans,sys},varargin);
    end

end

function aux_checkInputArgs(name,inv,trans,sys,n_in)

if CHECKS_ENABLED && n_in > 0

    inputArgsCheck({{name,'att',{'char','string'}};...
                    {inv,'att','contSet'};...
                    {trans,'att',{'transition','cell'}};...  % cell only legacy
                    {sys,'att','contDynamics'}});

end

end

function [name,inv,trans,sys] = aux_computeProperties(name,inv,trans,sys)

% loc.transition has to be an array of transition objects
if ~isa(trans,'transition')
    if iscell(trans)
        % legacy error message
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Transitions have to be an object array instead of a cell array.'));
    else
        throw(CORAerror('CORA:wrongInputInConstructor',...
            'Transitions have to be a transition object array.'));
    end
end

% convert invariant sets to polytopes if possible (keep polytopes,
% levelSet, and fullspace/emptySet invariants)
if ~(isa(inv,'fullspace') || isa(inv,'polytope') ...
        || isa(inv,'levelSet') || isa(inv,'emptySet'))
    inv = polytope(inv);
end

end

% ------------------------------ END OF CODE ------------------------------
