function child = InstantiateBC(parent,template,bind)
% InstantiateBC - Apply a bind to a Base-Component template, this
%    instantiates the component using the variable names from the bind
%    in the invariant/flow/guard/reset equations; additionally, all
%    constant parameters are removed from the list of variables
%
% Syntax:  
%    child = InstantiateBC(parent,template,bind)
%
% Inputs:
%    parent (struct) - instantiated parent component
%    template (struct) - template of BC to instantiate
%    bind (struct) - bind (i.e., instantion in parent network component)
%                    to be used for instantiation of the base component
%
% Outputs:
%    child (struct) - instantiated base component (child of "parent"
%                     network component)
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       ???
% Written:      ???
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% start with a copy of the template
child = template;

% save global name for new instance
child.name = parent.name + "." + bind.localName;

% output on command window
fprintf("instantiating Base Component: '%s'\n",child.name);

% Apply variable mapping to all equations

% number of locations in base components
numStates = length(child.States);

% loop over each location
for st = 1:numStates
    % read out location
    St = child.States(st);
    
    % Apply mapping to Flow
    child.States(st).Flow = applyMappingToEquation(St.Flow,bind);
        
    % Apply Mapping to Invariant
    child.States(st).Invariant = applyMappingToCondition(St.Invariant,bind);
    
    % number of outgoing transitions
    numTrans = length(St.Trans);

    % loop over each outgoing transition
    for tr = 1:numTrans
        % read out transition
        Tr = St.Trans(tr);
        
        % Apply mapping to reset function
        child.States(st).Trans(tr).reset = applyMappingToEquation(Tr.reset,bind);
        
        % Apply mapping to guard set
        child.States(st).Trans(tr).guard = applyMappingToCondition(Tr.guard,bind);
    end
end

% init array that all variables are states (and not constants)
noConst = true(numel(child.listOfVar),1);

% loop over each variable
for v=1:numel(child.listOfVar)
    
    % index of variable to be renamed
    renameIdx = bind.keys == child.listOfVar(v).name;

    if any(renameIdx)
        % rename state/constant
        child.listOfVar(v).name = bind.renames(renameIdx);

        % check if state is constant (converting a string that does not 
        % represent a double value to a double results in NaN)
        if ~isnan(str2double(bind.renames(renameIdx)))
            noConst(v) = false;
        end
    end

end

% remove all constants from list of variables
child.listOfVar = child.listOfVar(noConst);

%------------- END OF CODE --------------
