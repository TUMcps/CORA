classdef parallelHybridAutomaton
% hybridAutomaton - Object and Copy Constructor 
%
% Syntax:  
%    obj = parallelHybridAutomaton(components,stateBinds,inputBinds)
%
% Inputs:
%    components - cell array of hybridAutomaton objects that represent the
%                 subcomponents
%    stateBinds - cell array of nx1 int arrays. Maps component states to 
%                 states of composed system
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
    components = [];
    
    % description of the composition of the single subcomponents
    bindsStates = [];
    
    % description of the connection of the single subcomponents
    bindsInputs = [];
    
    % struct storing the results of simulation and reachability analysis
    result = [];  
    
    % number of states for the complete automaton
    numStates = [];
    
    % number of inputs for the complete automaton
    numInputs = [];
end

methods
    
    % Class Constructor
    function Obj = parallelHybridAutomaton(varargin)
        
        % no argument is passed (default constructor)
        if nargin == 0
            
            disp('ParallelHybridAutomaton needs more input values');
            Obj=[];


        elseif nargin == 3
            
            % parse input arguments
            Obj.components=varargin{1};

            stateBinds = varargin{2};
            inputBinds = varargin{3};

            % make sure cells are column arrays
            if size(stateBinds,2) > 1
                stateBinds = stateBinds.';
            end
            if size(inputBinds,2) > 1
                inputBinds = inputBinds.';
            end
            
            % extract overall state and input dimensions
            [Obj.numStates,Obj.numInputs] = validateBinds(stateBinds,inputBinds);
            
            % initialize remaining properties
            Obj.bindsStates = stateBinds;
            Obj.bindsInputs = inputBinds;

            Obj.result=[];


        % not enough or too many inputs are passed    
        else
            disp('This class needs more/less input values');
            Obj=[];
        end
    end
end
end

%------------- END OF CODE --------------