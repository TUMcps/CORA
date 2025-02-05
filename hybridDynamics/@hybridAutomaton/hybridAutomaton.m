classdef hybridAutomaton < hybridDynamics
% hybridAutomaton - constructor for class hybridAutomaton
%
% Syntax:
%    HA = hybridAutomaton()
%    HA = hybridAutomaton(loc)
%    HA = hybridAutomaton(name,loc)
%
% Inputs:
%    loc - location-array storing location objects
%
% Outputs:
%    name - name of automaton
%    HA - generated hybridAutomaton object
%
% Example:
%    % invariant
%    inv = polytope([-1,0],0);
% 
%    % transition
%    guard = polytope([0,1],0,[-1,0],0);
%    reset = linearReset([1,0;0,-0.75]);
%    trans(1) = transition(guard,reset,1);
% 
%    % flow equation
%    dynamics = linearSys([0,1;0,0],[0;0],[0;-9.81]);
%
%    % define location
%    loc(1) = location('S1',inv,trans,dynamics);
% 
%    % instantiate hybrid automaton
%    HA = hybridAutomaton('bouncingball',loc);
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
%                15-October-2024 (MW, add name property)
%                16-October-2024 (TL, renames dim to nrOfStates)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = private, GetAccess = public)
    name = '';                  % name of automaton
    location = location();	    % cell-array of location objects

    % internally-set properties
    nrOfDims;                   % number of states of each location
    nrOfInputs;                 % number of inputs for each location
    nrOfOutputs;                % number of outputs for each location
    nrOfDisturbances;           % number of disturbances for each location
    nrOfNoises;                 % number of noises for each location

    % legacy
    dim;                        % (also) state dimension
    nrOfStates;                 % (also) state dimension
end

methods
    
    % class constructor
    function HA = hybridAutomaton(varargin)

        % 0. empty
        assertNarginConstructor(0:2,nargin);
        if nargin == 0
            return
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'hybridAutomaton')
            HA = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [name,locs] = aux_parseInputArgs(varargin{:});

        % 3. check correctness of input arguments
        aux_checkInputArgs(name,locs,nargin);

        % 4. compute internal properties
        [states,inputs,outputs,dists,noises] = aux_computeProperties(locs);

        % 5. assign properties
        HA.name = name;
        HA.location = reshape(locs,[],1);
        HA.nrOfDims = states;
        HA.nrOfInputs = inputs;
        HA.nrOfOutputs = outputs;
        HA.nrOfDisturbances = dists;
        HA.nrOfNoises = noises;

    end
end


% getter & setter ---------------------------------------------------------

methods 
    function dim = get.dim(sys)
        CORAwarning('CORA:deprecated', 'property', 'contDynamics.dim', 'CORA v2025', ...
            'Please use contDynamics.nrOfDims instead.', ...
            'This change was made to be consistent with the other properties.')
        dim = sys.nrOfDims;
    end
    function sys = set.dim(sys,dim)
        CORAwarning('CORA:deprecated', 'property', 'contDynamics.dim', 'CORA v2025', ...
            'Please use contDynamics.nrOfDims instead.', ...
            'This change was made to be consistent with the other properties.')
        sys.nrOfDims = dim;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [name,locs] = aux_parseInputArgs(varargin)
% nargin may only be 1 or 2

name = 'hybridAutomaton';
if nargin == 1
    locs = varargin{1};
elseif nargin == 2
    [name,locs] = varargin{:};
end

end

function aux_checkInputArgs(name,locs,n_in)

    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({{name,'att',{'char','string'}};...
                        {locs,'att',{'location','cell'}}});  % cell checked below (legacy)...

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

            % 1. invariant of each location must have same dimension as
            % flow equation of that same location (unless empty)
            if dim(locs(i).invariant) ~= locs(i).contDynamics.nrOfDims
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Ambient dimension of invariant has to match state dimension of flow.'));
            end

            % 2. guard set of each transition of a location must have
            % same dimension as flow equation of that same location
            for j=1:length(locs(i).transition)
                if ~isnumeric(locs(i).transition(j).guard) && ...
                        dim(locs(i).transition(j).guard) ~= locs(i).contDynamics.nrOfDims
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Ambient dimension of guard set has to match state dimension of flow.'));
                end
            end

            % 3. target of each transition must be <= number of locations
            % (emptiness, integer value, and larger than 0 already checked)
            for j=1:length(locs(i).transition)
                if any(locs(i).transition(j).target > numLoc)
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        'Targets exceed number of locations.'));
                end
            end

            % 4. output dimension of reset function must have same
            % dimension as flow equation of target dimension
            for j=1:length(locs(i).transition)
                % skip empty transitions
                if ~isemptyobject(locs(i).transition(j)) ...
                        && locs(i).transition(j).reset.postStateDim ...
                        ~= locs(locs(i).transition(j).target).contDynamics.nrOfDims
                    throw(CORAerror('CORA:wrongInputInConstructor',...
                        ['Output dimension of reset function has to match '...
                         'the state dimension of the flow equation of the target location.']));
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

function [states,inputs,outputs,dists,noises] = aux_computeProperties(locs)
% compute the number of states, inputs, outputs, disturbances, and noises
% for each location

    % loop over flows of all locations
    states = arrayfun(@(x) x.contDynamics.nrOfDims,locs,'UniformOutput',true);
    inputs = arrayfun(@(x) x.contDynamics.nrOfInputs,locs,'UniformOutput',true);
    outputs = arrayfun(@(x) x.contDynamics.nrOfOutputs,locs,'UniformOutput',true);
    dists = arrayfun(@(x) x.contDynamics.nrOfDisturbances,locs,'UniformOutput',true);
    noises = arrayfun(@(x) x.contDynamics.nrOfNoises,locs,'UniformOutput',true);

end

% ------------------------------ END OF CODE ------------------------------
