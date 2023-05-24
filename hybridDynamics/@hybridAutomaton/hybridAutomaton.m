classdef hybridAutomaton
% hybridAutomaton - constructor for class hybridAutomaton
%
% Syntax:
%    HA = hybridAutomaton()
%    HA = hybridAutomaton(loc)
%
% Inputs:
%    loc - location-array storing location objects
%
% Outputs:
%    HA - generated hybridAutomaton object
%
% Example:
%    % invariant
%    inv = mptPolytope([-1,0],0);
% 
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%    trans = transition(guard,reset,1);
% 
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
%
%    % define location
%    loc = location('S1',inv,trans,dynamics);
% 
%    % instantiate hybrid automaton
%    HA = hybridAutomaton(loc);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, location, transition

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      03-May-2007 
% Last update:  16-June-2022 (MW, add checks for object properties)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = private, GetAccess = public)
    % cell-array of location objects
    location = location();
end

methods
    
    % class constructor
    function HA = hybridAutomaton(varargin)

        % parse input arguments
        if nargin == 0
            return
        elseif nargin == 1 && isa(varargin{1},'hybridAutomaton')
            % copy constructor
            HA = varargin{1}; return
        elseif nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs',1));
        elseif ~isa(varargin{1},'location')
            throw(CORAerror('CORA:wrongInputInConstructor',...
                'Input argument has to be an array of location objects.'));
        end
        
        HA.location = varargin{1};
        
        % check input arguments
        if CHECKS_ENABLED

            % number of locations
            numLoc = length(HA.location);
            
            for i=1:numLoc
                % 1. invariant of each location has to have same dimension as
                % flow equation of that same location (unless empty)
                if ~(isnumeric(HA.location(i).invariant) && ...
                        isempty(HA.location(i).invariant))
                    if dim(HA.location(i).invariant) ~= HA.location(i).contDynamics.dim
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Ambient dimension of invariant has to match state dimension of flow.'));
                    end
                end
    
                % 2. guard set of each transition of a location has to have
                % same dimension as flow equation of that same location
                for j=1:length(HA.location(i).transition)
                    if ~isnumeric(HA.location(i).transition(j).guard) && ...
                            dim(HA.location(i).transition(j).guard) ~= HA.location(i).contDynamics.dim
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Ambient dimension of guard set has to match state dimension of flow.'));
                    end
                end
    
                % 3. target of each transition has to be <= number of locations
                % (emptiness, integer value, and larger than 0 already checked)
                for j=1:length(HA.location(i).transition)
                    if any(HA.location(i).transition(j).target > numLoc)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            'Targets exceed number of locations.'));
                    end
                end
    
                % 4. output dimension of reset function has to have same
                % dimension as flow equation of target dimension
                for j=1:length(HA.location(i).transition)
                    msg = ['Output dimension of reset function has to match '...
                            'the state dimension of the flow equation of the target location.'];
                    if isfield(HA.location(i).transition(j).reset,'A') && ...
                            size(HA.location(i).transition(j).reset.A,1) ...
                                ~= HA.location(HA.location(i).transition(j).target).contDynamics.dim
                        % linear reset function -> A checked
                        throw(CORAerror('CORA:wrongInputInConstructor',msg));
                    elseif isfield(HA.location(i).transition(j).reset,'f')
                        % nonlinear reset function -> output dim of f checked
                        [~,nrOutputStates] = inputArgsLength(HA.location(i).transition(j).reset.f);
                        if nrOutputStates ~= HA.location(HA.location(i).transition(j).target).contDynamics.dim
                            throw(CORAerror('CORA:wrongInputInConstructor',msg));
                        end
                    end
                end
    
                % 5. no duplicates in synchonization labels of a location
                % (excluded case where no synchronization label is given)
                syncLabelList = {};
                for j=1:length(HA.location(i).transition)
                    syncLabel = HA.location(i).transition(j).syncLabel;
                    if ~isempty(syncLabel)
                        if ismember(syncLabel,syncLabelList)
                            throw(CORAerror('CORA:wrongInputInConstructor',...
                                ['Each synchonization label may only be used ...' ...
                                'in one outgoing transition per location.']));
                        end
                        % extend list of synchronization labels
                        syncLabelList = [syncLabelList;...
                            HA.location(i).transition(j).syncLabel];
                    end
                end
    
                % 6. unless synchronization labels differ, no more than one
                % immediate transition per location allowed
                emptyGuardSets = 0;
                for j=1:length(HA.location(i).transition)
                    if isnumeric(HA.location(i).transition(j).guard) && ...
                            isempty(HA.location(i).transition(j).guard) && ...
                            isempty(HA.location(i).transition(j).syncLabel)
                        % immediate transition without synchronization label
                        % (since duplicates of synchronization labels already
                        % checked before, it suffices to check cases where no
                        % synchronization label is given)
                        emptyGuardSets = emptyGuardSets + 1;
                    end
                    if emptyGuardSets > 1
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            ['Once one immediate transition (guard = []) per location allowed,'...
                             'unless synchronization labels differ.']));
                    end
                end
    
            end
        end

    end
end
end

%------------- END OF CODE --------------