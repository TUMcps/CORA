classdef parallelHybridAutomaton
% parallelHybridAutomaton - Object and Copy Constructor 
%
% Syntax:  
%    pHA = parallelHybridAutomaton(components,inputBinds)
%
% Inputs:
%    components - cell array of hybridAutomaton objects that represent the
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
%    inv = mptPolytope(1,T_off);
%    guard = conHyperplane(1,T_off);
%    reset = struct('A',1,'c',0);
%    trans = {transition(guard,reset,2)};
%    linSys = linearSys(-(a1 + b1),[b1,a1],c1,1);
%    loc{1} = location('on',inv,trans,linSys);
%
%    % first component, second location
%    inv = mptPolytope(-1,-T_on);
%    guard = conHyperplane(1,T_on);
%    reset = struct('A',1,'c',0);
%    trans = {transition(guard,reset,1)};
%    linSys = linearSys(-(a1 + b1),[b1,a1],[],1);
%    loc{2} = location('off',inv,trans,linSys);
%
%    HA1 = hybridAutomaton(loc);
%
%    % second component, first location
%    inv = mptPolytope(1,T_off);
%    guard = conHyperplane(1,T_off);
%    reset = struct('A',1,'c',0);
%    trans = {transition(guard,reset,2)};
%    linSys = linearSys(-(a2 + b2),[b2,a2],c2,1);
%    loc{1} = location('on',inv,trans,linSys);
%
%    % second component, second location
%    inv = mptPolytope(-1,-T_on);
%    guard = conHyperplane(1,T_on);
%    reset = struct('A',1,'c',0);
%    trans = {transition(guard,reset,1)};
%    linSys = linearSys(-(a2 + b2),[b2,a2],[],1);
%    loc{2} = location('off',inv,trans,linSys);
%
%    HA2 = hybridAutomaton(loc);
%
%    % parallel hybrid automaton
%    components = {HA1,HA2};
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
% See also: ---

% Author:       Johann Schoepfer, Niklas Kochdumper, Mark Wetzlinger
% Written:      05-June-2018
% Last update:  16-June-2022 (MW, add input argument checks)
% Last revision:---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    
    % list of hybridAutomaton objects representing the subcomponents of the
    % system
    components (:,1) cell = {};
    
    % description of the connection of the single subcomponents
    bindsInputs (:,1) cell = {};

    % --- the remaining properties are set internally ---
    % description of the composition of the single subcomponents
    bindsStates (:,1) cell = {};
    
    % number of states for the complete automaton
    numStates (1,1) {mustBeNonnegative,mustBeFinite} = 0;
    
    % number of inputs for the complete automaton
    numInputs (1,1) {mustBeNonnegative,mustBeFinite} = 0;
end

methods
    
    % Class Constructor
    function pHA = parallelHybridAutomaton(varargin)

        % throw error for wrong number of inputs
        if nargin < 2
            throw(CORAerror('CORA:notEnoughInputArgs',2));
        elseif nargin > 2
            throw(CORAerror('CORA:tooManyInputArgs',2));
        end
            
        % assign subcomponents of parallel hybrid automaton
        pHA.components = varargin{1};
        numComps = length(pHA.components);
        
        % check whether each component is a hybrid automaton
        if ~all(cellfun(@(x) isa(x,'hybridAutomaton'),pHA.components,...
                'UniformOutput',true))
            throw(CORAerror('CORA:wrongInputInConstructor',...
                ['Each component of the parallel hybrid automaton '...
                'has to be a hybridAutomaton object.']));
        end

        % check locations of each component
        for i=1:numComps
            % each location has to have same number of states
            states = cellfun(@(x) length(x.contDynamics.dim),...
                pHA.components{i}.location,'UniformOutput',true);
            if ~all(states == states(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of states.']));
            end

            % each location has to have same number of inputs
            inputs = cellfun(@(x) length(x.contDynamics.nrOfInputs),...
                pHA.components{1}.location,'UniformOutput',true);
            if ~all(inputs == inputs(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of inputs.']));
            end

            % each location has to have same number of outputs
            outputs = cellfun(@(x) length(x.contDynamics.nrOfOutputs),...
                pHA.components{1}.location,'UniformOutput',true);
            if ~all(outputs == outputs(1))
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Within a component of the parallel hybrid automaton, '...
                    'the continuous dynamics of each location have to '...
                    'have the same number of outputs.']));
            end
        end

        % assign input binds
        if size(varargin{2},2) > 1
            pHA.bindsInputs = varargin{2}.';
        else
            pHA.bindsInputs = varargin{2};
        end

        % check input binds (nx2 array, integer, ...)
        for i=1:length(pHA.bindsInputs)
            % nrOfInputs(max)-by-2 array
            % (we can use any location since all must have same nrOfInputs)
            if ( isempty(pHA.bindsInputs{i}) ...
                    && pHA.components{i}.location{1}.contDynamics.nrOfInputs ~= 0 ) ...
                    || ( size(pHA.bindsInputs{i},2) ~= 2 || ...
                    size(pHA.bindsInputs{i},1) > pHA.components{i}.location{1}.contDynamics.nrOfInputs )
                throw(CORAerror('CORA:wrongInputInConstructor',...
                    ['Input binds have to be\n'...
                     '   empty if that component''s number of inputs is zero\n'...
                     '   an mx2 array, where m is the number of inputs of that component.']));
            end
            % integer
            inputBindsArray = pHA.bindsInputs{i};
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

        % extract overall state and input dimensions (remaining properties)
        pHA = validateBinds(pHA);
    end
end
end

%------------- END OF CODE --------------