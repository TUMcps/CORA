classdef parallelHybridAutomaton
% hybridAutomaton - Object and Copy Constructor 
%
% Syntax:  
%    obj = parallelHybridAutomaton(components,inputBinds)
%
% Inputs:
%    components - cell array of hybridAutomaton objects that represent the
%                 subcomponents
%    inputBinds - cell array of nx2 int arrays. Maps component inputs to 
%                 states/inputs of composed system
%
% Outputs:
%    obj - generated parallelHybridAutomaton object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author: Johann Schoepfer, Niklas Kochdumper
% Written: 05-June-2018
% Last update: ---
% Last revision: ---

%------------- BEGIN CODE --------------

properties (SetAccess = protected, GetAccess = public)
    
    % list of hybridAutomaton objects representing the subcomponents of the
    % system
    components (:,1) {cell} = {}
    
    % description of the composition of the single subcomponents
    bindsStates (:,1) {cell} = {}
    
    % description of the connection of the single subcomponents
    bindsInputs (:,1) {cell} = {}
    
    % number of states for the complete automaton
    numStates (1,1) {mustBeNonnegative,mustBeFinite} = 0;
    
    % number of inputs for the complete automaton
    numInputs (1,1) {mustBeNonnegative,mustBeFinite} = 0;
end

methods
    
    % Class Constructor
    function Obj = parallelHybridAutomaton(components,inputBinds)
            
        % parse input arguments
        Obj.components = components;
        
        if size(inputBinds,2) > 1
            Obj.bindsInputs = inputBinds.';
        else
            Obj.bindsInputs = inputBinds;
        end
        
        for i = 1:length(Obj.components)
           if ~isa(Obj.components{i},'hybridAutomaton')
              error('Wrong format for input argument "components"!'); 
           end
        end

        % extract overall state and input dimensions
        Obj = validateBinds(Obj);
    end
end
end

%------------- END OF CODE --------------