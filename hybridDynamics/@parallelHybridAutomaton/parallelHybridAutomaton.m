classdef parallelHybridAutomaton < hybridDynamics
% parallelHybridAutomaton - Object and Copy Constructor 
%
% Syntax:
%    pHA = parallelHybridAutomaton()
%    pHA = parallelHybridAutomaton(components,inputBinds)
%    pHA = parallelHybridAutomaton(name,components,inputBinds)
%
% Inputs:
%    name - name of the parallel hybrid automaton
%    components - class array of hybridAutomaton objects representing the
%                 subcomponents
%    inputBinds - cell array of nx2 int arrays. Maps component inputs to 
%                 states/inputs of composed system
%
% Outputs:
%    pHA - generated parallelHybridAutomaton object
%
% Example:
%    % parameters
%    a1 = 0.5; b1 = 0.4; c1 = 6;
%    a2 = 0.5; b2 = 0.3; c2 = 7;
%    T_off = 21; T_on = 20;
%
%    % first component, first location
%    inv = polytope(1,T_off);
%    guard = polytope([],[],1,T_off);
%    reset = linearReset(1);
%    trans(1) = transition(guard,reset,2);
%    linsys = linearSys(-(a1 + b1),[b1,a1],c1,1);
%    loc(1) = location('on',inv,trans,linsys);
%
%    % first component, second location
%    inv = polytope(-1,-T_on);
%    guard = polytope([],[],1,T_on);
%    reset = linearReset(1);
%    trans(1) = transition(guard,reset,1);
%    linsys = linearSys(-(a1 + b1),[b1,a1],[],1);
%    loc(2) = location('off',inv,trans,linsys);
%
%    HA1 = hybridAutomaton(loc);
%
%    % second component, first location
%    inv = polytope(1,T_off);
%    guard = polytope([],[],1,T_off);
%    reset = linearReset(1);
%    trans(1) = transition(guard,reset,2);
%    linsys = linearSys(-(a2 + b2),[b2,a2],c2,1);
%    loc(1) = location('on',inv,trans,linsys);
%
%    % second component, second location
%    inv = polytope(-1,-T_on);
%    guard = polytope([],[],1,T_on);
%    reset = linearReset(1);
%    trans(1) = transition(guard,reset,1);
%    linsys = linearSys(-(a2 + b2),[b2,a2],[],1);
%    loc(2) = location('off',inv,trans,linsys);
%
%    HA2 = hybridAutomaton(loc);
%
%    % parallel hybrid automaton
%    components = [HA1,HA2];
%    inputBinds{1} = [0 1; ...   % first global input
%                     2 1];      % first output of component 2
%    inputBinds{2} = [0 1; ...   % first global input
%                     1 1];      % first output of component 1
%    pHA = parallelHybridAutomaton('thermostat',components,inputBinds);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: hybridAutomaton, location

% Authors:       Johann Schoepfer, Niklas Kochdumper, Mark Wetzlinger
% Written:       05-June-2018
% Last update:   16-June-2022 (MW, add input argument checks)
%                06-June-2023 (MW, add mergedTrans property)
%                21-June-2023 (MW, add/rename dim/nrOfInputs/nrOfOutputs)
%                15-October-2024 (MW, add name property)
%                16-October-2024 (TL, renames dim to nrOfStates)
% Last revision: 21-June-2023 (MW, restructure, integrate validateBinds)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    
    % --- properties provided by user ---
    name = '';                          % name of the automaton
    components = hybridAutomaton();     % subcomponents of the system
    bindsInputs = {};                   % connection of subcomponents

    % --- properties computed from user inputs ---
    bindsStates = {};                   % mapping local -> global state
    nrOfStates;                         % number of states of composed automaton
    nrOfInputs;                         % number of global inputs
    nrOfOutputs;                        % number of global outputs of composed automaton
    nrOfDisturbances;                   % number of global disturbances of composed automaton
    nrOfNoises;                         % number of global noises of composed automaton

    % --- properties stored on-the-fly for faster evaluation ---
    locProd;                            % computed location products
    mergedTrans;                        % merged transitions

    % legacy
    dim;                % (also) state dimension
    
end

methods
    
    % Class Constructor
    function pHA = parallelHybridAutomaton(varargin)

        % 0. empty
        assertNarginConstructor(0:3,nargin);
        if nargin == 0
            return;
        end

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'parallelHybridAutomaton')
            pHA = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        [name,comps,inputBinds] = aux_parseInputArgs(varargin{:});
            
        % 3. check correctness of input arguments
        aux_checkInputArgs(name,comps,inputBinds,nargin);
            
        % 4. compute internal properties
        [name,comps,inputBinds,stateBinds,states,inputs,outputs,dists,noises] = ...
            aux_computeProperties(name,comps,inputBinds);

        % 5. assign properties
        pHA.name = name;
        pHA.components = comps;
        pHA.bindsInputs = inputBinds;
        pHA.bindsStates = stateBinds;
        pHA.nrOfStates = states;
        pHA.nrOfInputs = inputs;
        pHA.nrOfOutputs = outputs;
        pHA.nrOfDisturbances = dists;
        pHA.nrOfNoises = noises;

        % init structs for computed location products/merged transitions
        pHA.locProd = struct('location',cell(0),'locID',cell(0));
        pHA.mergedTrans = struct('transition',cell(0),'locID',cell(0),'transID',cell(0));
    end

end


% getter & setter ---------------------------------------------------------

methods 
    function dim = get.dim(sys)
        CORAwarning('CORA:deprecated', 'property', 'contDynamics.dim', 'CORA v2025', ...
            'Please use contDynamics.nrOfStates instead.', ...
            'This change was made to be consistent with the other properties.')
        dim = sys.nrOfStates;
    end
    function sys = set.dim(sys,dim)
        CORAwarning('CORA:deprecated', 'property', 'contDynamics.dim', 'CORA v2025', ...
            'Please use contDynamics.nrOfStates instead.', ...
            'This change was made to be consistent with the other properties.')
        sys.nrOfStates = dim;
    end
end

end


% Auxiliary functions -----------------------------------------------------

function [name,comps,inputBinds] = aux_parseInputArgs(varargin)
    
    narginchk(2,3);

    % default name
    name = 'parallelHybridAutomaton';
    if nargin == 2
        [comps,inputBinds] = varargin{:};
    elseif nargin == 3
        [name,comps,inputBinds] = varargin{:};
    end

end

function aux_checkInputArgs(name,comps,inputBinds,n_in)

    if CHECKS_ENABLED && n_in > 0

        inputArgsCheck({{name,'att',{'char','string'}};...
                        {comps,'att','hybridAutomaton','vector'};...
                        {inputBinds,'att','cell'}});

        % read out number of components
        numComps = numel(comps);
        
        % check whether each component is a hybrid automaton
        if ~all(arrayfun(@(x) isa(x,'hybridAutomaton'),comps,...
                'UniformOutput',true))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Each component of the parallel hybrid automaton '...
                'has to be a hybridAutomaton object.']));
        end

        % check locations of each component
        for i=1:numComps
            % each location has to have same number of states
            states = arrayfun(@(x) length(x.contDynamics.nrOfStates),...
                comps(i).location,'UniformOutput',true);
            if ~all(states == states(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of states.']));
            end

            % each location has to have same number of inputs
            inputs = arrayfun(@(x) length(x.contDynamics.nrOfInputs),...
                comps(i).location,'UniformOutput',true);
            if ~all(inputs == inputs(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of inputs.']));
            end

            % each location has to have same number of outputs
            outputs = arrayfun(@(x) length(x.contDynamics.nrOfOutputs),...
                comps(i).location,'UniformOutput',true);
            if ~all(outputs == outputs(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of outputs.']));
            end

            % each location has to have the same number of disturbances
            dists = arrayfun(@(x) length(x.contDynamics.nrOfDisturbances),...
                comps(i).location,'UniformOutput',true);
            if ~all(dists == dists(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of disturbances.']));
            end

            % each location has to have the same number of noises
            noises = arrayfun(@(x) length(x.contDynamics.nrOfNoises),...
                comps(i).location,'UniformOutput',true);
            if ~all(noises == noises(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of noises.']));
            end
        end

        % check input binds (nx2 array, integer, ...)
        for i=1:numel(inputBinds)
            % nrOfInputs(max)-by-2 array
            % (we can use any location since all must have same nrOfInputs)
            if xor(isempty(inputBinds{i}), ...
                    comps(i).location(1).contDynamics.nrOfInputs == 0) ...
                    || ~any(size(inputBinds{i},2) ~= [0,2]) ...
                    || size(inputBinds{i},1) ~= comps(i).location(1).contDynamics.nrOfInputs
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Input binds have to be\n'...
                     '   empty if that component''s number of inputs is zero\n'...
                     '   an mx2 array, where m is the number of inputs of that component.']));
            end
            % integer
            inputBindsArray = inputBinds{i};
            if any(any(mod(inputBindsArray,1)))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Input binds have to consist of integer values.'));
            end
            % left entry: global or any other component
            if any(inputBindsArray(:,1) < 0) || any(inputBindsArray(:,1) > numComps)
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    'Left column of input bind has to be >= 0 and <= number of components.'));
            end
            
        end
    end

end

function [name,comps,inputBinds,stateBinds,states,inputs,outputs,dists,noises] = ...
    aux_computeProperties(name,comps,inputBinds)
% determines the internal properties stateBinds, number of states/global
% inputs/global outputs of the composed automaton
% note: parts were previously handled in validateBinds (now removed)

    % reshape to column vectors
    comps = reshape(comps,[],1);
    inputBinds = reshape(inputBinds,[],1);

    % determine number of states of composed automaton
    n_comp = arrayfun(@(x) x.location(1).contDynamics.nrOfStates,comps,'UniformOutput',true)';
    states = sum(n_comp);

    % compute state binds
    stateBinds = cell(length(comps),1);
    n_end = cumsum(n_comp);
    n_start = [1; reshape(n_end(1:end-1)+1,[],1)];
    for i=1:length(stateBinds)
        if n_start(i) ~= n_end(i)
            stateBinds{i} = (n_start(i):n_end(i))';
        else
            stateBinds{i} = n_start(i);
        end
    end

    % determine number of global inputs to composed automaton
    % concatenate bind arrays vertically
    concatInputs = vertcat(inputBinds{:});
    
    % check whether inputBinds are valid
    inputComps = concatInputs(:,1);
    inputIndices = concatInputs(:,2);
    
    % get number of global inputs
    uniqueIndices = unique(inputIndices(inputComps == 0));
    inputs = length(uniqueIndices);

    % determine number of global outputs to composed automaton: since CORA
    % does not support output equations for composed parallel hybrid
    % automata, we assume an identity output matrix, i.e., y = x
    outputs = states;

    % determine number of global disturbances and global noises:
    % both disturbances and noises per component are independent from one
    % another; also, they are by input argument check forced to be the same
    % for each location, so we can just use the first one for each comp
    dists = sum(arrayfun(@(x) x.location(1).contDynamics.nrOfDisturbances, comps));
    noises = sum(arrayfun(@(x) x.location(1).contDynamics.nrOfNoises, comps));

end

% ------------------------------ END OF CODE ------------------------------
