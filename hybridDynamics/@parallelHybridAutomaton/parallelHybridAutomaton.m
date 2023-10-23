classdef parallelHybridAutomaton < hybridDynamics
% parallelHybridAutomaton - Object and Copy Constructor 
%
% Syntax:
%    pHA = parallelHybridAutomaton()
%    pHA = parallelHybridAutomaton(components,inputBinds)
%
% Inputs:
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
%    guard = conHyperplane(1,T_off);
%    reset = struct('A',1,'c',0);
%    trans(1) = transition(guard,reset,2);
%    linSys = linearSys(-(a1 + b1),[b1,a1],c1,1);
%    loc(1) = location('on',inv,trans,linSys);
%
%    % first component, second location
%    inv = polytope(-1,-T_on);
%    guard = conHyperplane(1,T_on);
%    reset = struct('A',1,'c',0);
%    trans(1) = transition(guard,reset,1);
%    linSys = linearSys(-(a1 + b1),[b1,a1],[],1);
%    loc(2) = location('off',inv,trans,linSys);
%
%    HA1 = hybridAutomaton(loc);
%
%    % second component, first location
%    inv = polytope(1,T_off);
%    guard = conHyperplane(1,T_off);
%    reset = struct('A',1,'c',0);
%    trans(1) = transition(guard,reset,2);
%    linSys = linearSys(-(a2 + b2),[b2,a2],c2,1);
%    loc(1) = location('on',inv,trans,linSys);
%
%    % second component, second location
%    inv = polytope(-1,-T_on);
%    guard = conHyperplane(1,T_on);
%    reset = struct('A',1,'c',0);
%    trans(1) = transition(guard,reset,1);
%    linSys = linearSys(-(a2 + b2),[b2,a2],[],1);
%    loc(2) = location('off',inv,trans,linSys);
%
%    HA2 = hybridAutomaton(loc);
%
%    % parallel hybrid automaton
%    components = [HA1,HA2];
%    inputBinds{1} = [0 1; ...   % first global input
%                     2 1];      % first output of component 2
%    inputBinds{2} = [0 1; ...   % first global input
%                     1 1];      % first output of component 1
%    pHA = parallelHybridAutomaton(components,inputBinds);
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
% Last revision: 21-June-2023 (MW, restructure, integrate validateBinds)

% ------------------------------ BEGIN CODE -------------------------------

properties (SetAccess = protected, GetAccess = public)
    
    % list of hybridAutomaton objects representing the subcomponents of the
    % system
    components = hybridAutomaton();
    
    % description of the connection of the single subcomponents
    bindsInputs = {};

    % --- the remaining properties are set internally ---
    % description of the composition of the single subcomponents
    bindsStates = {};

    % number of states of composed automaton
    dim;

    % number of global inputs
    nrOfInputs;

    % number of global outputs of composed automaton
    nrOfOutputs;

    % computed location products (used for speed up)
    locProd;

    % merged transitions (used for speed up)
    mergedTrans;
end

methods
    
    % Class Constructor
    function pHA = parallelHybridAutomaton(varargin)

        % 1. copy constructor
        if nargin == 1 && isa(varargin{1},'parallelHybridAutomaton')
            pHA = varargin{1}; return
        end

        % 2. parse input arguments: varargin -> vars
        if nargin == 0
            return;
        elseif nargin < 2
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs',2));
        else
            comps = varargin{1};
            % assign input binds
            if size(varargin{2},2) > 1
                inputBinds = varargin{2}.';
            else
                inputBinds = varargin{2};
            end
        end
            
        % 3. check correctness of input arguments
        aux_checkInputArgs(comps,inputBinds,nargin);
            
        % 4. compute internal properties
        [stateBinds,n,m,r] = aux_computeProperties(comps,inputBinds);

        % 5. assign properties
        pHA.components = comps;
        pHA.bindsInputs = inputBinds;
        pHA.bindsStates = stateBinds;
        pHA.dim = n;
        pHA.nrOfInputs = m;
        pHA.nrOfOutputs = r;

        % init structs for computed location products/merged transitions
        pHA.locProd = struct('location',cell(0),'locID',cell(0));
        pHA.mergedTrans = struct('transition',cell(0),'locID',cell(0),'transID',cell(0));
    end
end
end


% Auxiliary functions -----------------------------------------------------

function aux_checkInputArgs(comps,inputBinds,n_in)

    if CHECKS_ENABLED && n_in > 0

        % read out number of components
        numComps = length(comps);
        
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
            states = arrayfun(@(x) length(x.contDynamics.dim),...
                comps(i).location,'UniformOutput',true);
            if ~all(states == states(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of states.']));
            end

            % each location has to have same number of inputs
            inputs = arrayfun(@(x) length(x.contDynamics.nrOfInputs),...
                comps(1).location,'UniformOutput',true);
            if ~all(inputs == inputs(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of inputs.']));
            end

            % each location has to have same number of outputs
            outputs = arrayfun(@(x) length(x.contDynamics.nrOfOutputs),...
                comps(1).location,'UniformOutput',true);
            if ~all(outputs == outputs(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of outputs.']));
            end
        end

        % check input binds (nx2 array, integer, ...)
        for i=1:length(inputBinds)
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

function [stateBinds,n,m,r] = aux_computeProperties(comps,inputBinds)
% determines the internal properties stateBinds, number of states/global
% inputs/global outputs of the composed automaton
% note: parts were previously handled in validateBinds (now removed)

    % determine number of states of composed automaton
    n_comp = arrayfun(@(x) x.location(1).contDynamics.dim,comps,'UniformOutput',true)';
    n = sum(n_comp);

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
    temp = unique(inputIndices(inputComps == 0));
    m = length(temp);

    % determine number of global outputs to composed automaton: since CORA
    % does not support output equations for composed parallel hybrid
    % automata, we use the state dimension instead
    r = n;

end

% ------------------------------ END OF CODE ------------------------------
