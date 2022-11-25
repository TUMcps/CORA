function [newObj] = pruneDynamics(obj,options)
% pruneDynamics - removes part of the system dynamics, which is irrelevant
% for the reachable set of output values (this is only a prototype)
%
% Syntax:  
%    [newObj] = pruneDynamics(obj,options)
%
% Inputs:
%    obj - continuous system object
%    options - options for the computation of reachable sets
%
% Outputs:
%    newObj - pruned continuous system object
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Matthias Althoff
% Written:      04-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% instantiate newObj
newObj = obj;

% obtain incidence matrices
if isa(obj,'linearSys')
    % input to state incidence matrix
    inputToState = sparse(ones(obj.dim,obj.nrOfInputs));
    ind = obj.B==0;
    inputToState(ind) = 0;
    
    % state to state incidence matrix
    stateToState = sparse(ones(obj.dim,obj.dim));
    ind = obj.A==0;
    stateToState(ind) = 0;
    
    % input to output incidence matrix
    inputToOutput = sparse(ones(obj.nrOfOutputs,obj.nrOfInputs));
    ind = obj.D==0;
    inputToOutput(ind) = 0;    
    
    % state to output incidence matrix
    stateToOutput = sparse(ones(obj.nrOfOutputs,obj.dim));
    ind = obj.C==0;
    stateToOutput(ind) = 0; 
    
    % input incidence vector
    
    % initial state incidence vector
end


%% contributions
% state contributing to outputs
statesForOutputs = sum(stateToOutput)~=0;

% states influencing states that are relevant for the output
relevantStates = sum(stateToState(statesForOutputs,:))~=0;


%------------- END OF CODE --------------