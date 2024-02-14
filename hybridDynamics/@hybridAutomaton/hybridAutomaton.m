classdef hybridAutomaton < hybridDynamics
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
%    inv = polytope([-1,0],0);
% 
%    % transition
%    c = [-1;0]; d = 0; C = [0,1]; D = 0;
%    guard = conHyperplane(c,d,C,D);
%    reset = struct('A',[1,0;0,-0.75],'c',[0;0]);
%    trans(1) = transition(guard,reset,1);
% 
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
%
%    % define location
%    loc(1) = location('S1',inv,trans,dynamics);
% 
%    % instantiate hybrid automaton
%    HA = hybridAutomaton(loc);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: parallelHybridAutomaton, location, transition

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       03-May-2007 
% Last update:   16-June-2022 (MW, add checks for object properties)
%                21-June-2023 (MW, add internal properties, restructure)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    location = location();	    % cell-array of location objects

    % internally-set properties
    dim;                        % state dimension of each location
    nrOfInputs;                 % number of inputs for each location
    nrOfOutputs;                % number of outputs for each location
end

methods
    
    % class constructor
    function HA = hybridAutomaton(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'hybridAutomaton')
            HA = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        if nargin == 0
            return
        elseif nargin > 1
            throw(CORAerror('CORA:tooManyInputArgs',1));
        else
            locs = varargin{1};
        end

        % 3. check correctness of input arguments
        aux_checkInputArgs(locs,nargin);

        % 4. compute internal properties
        [n,m,r] = aux_computeProperties(locs);

        % 5. assign properties
        HA.location = locs;
        HA.dim = n;
        HA.nrOfInputs = m;
        HA.nrOfOutputs = r;        

    end
end

end


% Auxiliary functions -----------------------------------------------------

function aux_checkInputArgs(locs,n_in)

    if CHECKS_ENABLED && n_in > 0

        if ~isa(locs,'location')
            if iscell(locs)
                % legacy error message
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Locations have to be an object array instead of a cell array.'));
            else
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Locations have to be a location object array.'));
            end
        end

        % number of locations
        numLoc = length(locs);
        
        for i=1:numLoc

            % 1. invariant of each location has to have same dimension as
            % flow equation of that same location (unless empty)
            if dim(locs(i).invariant) ~= locs(i).contDynamics.dim
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Ambient dimension of invariant has to match state dimension of flow.'));
            end

            % 2. guard set of each transition of a location has to have
            % same dimension as flow equation of that same location
            for j=1:length(locs(i).transition)
                if ~isnumeric(locs(i).transition(j).guard) && ...
                        dim(locs(i).transition(j).guard) ~= locs(i).contDynamics.dim
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Ambient dimension of guard set has to match state dimension of flow.'));
                end
            end

            % 3. target of each transition has to be <= number of locations
            % (emptiness, integer value, and larger than 0 already checked)
            for j=1:length(locs(i).transition)
                if any(locs(i).transition(j).target > numLoc)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Targets exceed number of locations.'));
                end
            end

            % 4. output dimension of reset function has to have same
            % dimension as flow equation of target dimension
            for j=1:length(locs(i).transition)
                msg = ['Output dimension of reset function has to match '...
                        'the state dimension of the flow equation of the target location.'];
                if isfield(locs(i).transition(j).reset,'A') && ...
                        size(locs(i).transition(j).reset.A,1) ...
                            ~= locs(locs(i).transition(j).target).contDynamics.dim
                    % linear reset function -> A checked
                    throw(CORAerror('CORA:wrongInputInConstructor',msg));
                end
                if isfield(locs(i).transition(j).reset,'c') && ...
                        ~all(size(locs(i).transition(j).reset.c) ...
                            == [locs(locs(i).transition(j).target).contDynamics.dim,1])
                    % linear reset function -> dimension of c checked
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'The offset must be a column vector of proper dimension.'));
                end
                if isfield(locs(i).transition(j).reset,'f')
                    % nonlinear reset function -> output dim of f checked
                    [~,nrOutputStates] = inputArgsLength(locs(i).transition(j).reset.f);
                    if nrOutputStates ~= locs(locs(i).transition(j).target).contDynamics.dim
                        throw(CORAerror('CORA:wrongInputInConstructor',msg));
                    end
                end
            end
    
            % 5. no duplicates in synchonization labels of a location
            % (excluded case where no synchronization label is given)
            syncLabelList = {};
            for j=1:length(locs(i).transition)
                syncLabel = locs(i).transition(j).syncLabel;
                if ~isempty(syncLabel)
                    if ismember(syncLabel,syncLabelList)
                        throw(CORAerror('CORA:wrongInputInConstructor',...
                            ['Each synchonization label may only be used ...' ...
                            'in one outgoing transition per location.']));
                    end
                    % extend list of synchronization labels
                    syncLabelList = [syncLabelList;...
                        locs(i).transition(j).syncLabel];
                end
            end
    
            % 6. unless synchronization labels differ, no more than one
            % instant transition per location allowed
            emptyGuardSets = 0;
            for j=1:length(locs(i).transition)
                if isa(locs(i).transition(j).guard,'fullspace') && ...
                    isempty(locs(i).transition(j).syncLabel)
                    % instant transition without synchronization label
                    % (since duplicates of synchronization labels already
                    % checked before, it suffices to check cases where no
                    % synchronization label is given)
                    emptyGuardSets = emptyGuardSets + 1;
                end
                if emptyGuardSets > 1
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['Only one instant transition (guard = []) per location allowed,'...
                         'unless synchronization labels differ.']));
                end
    
            end
        end

    end

end

function [n,m,r] = aux_computeProperties(locs)
% compute the number of states, inputs, and outputs for each location

    % loop over flows of all locations
    n = arrayfun(@(x) x.contDynamics.dim,locs,'UniformOutput',true);
    m = arrayfun(@(x) x.contDynamics.nrOfInputs,locs,'UniformOutput',true);
    r = arrayfun(@(x) x.contDynamics.nrOfOutputs,locs,'UniformOutput',true);

end

% ------------------------------ END OF CODE ------------------------------
