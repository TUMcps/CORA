function [child] = InstantiateBC(parent,template,bind)
% Apply a bind to a Base-Component template
% INPUTS (all structs):
%   parent:     instantiated parent component
%   template:   template of BC to instantiate
%   bind:       bind to be used for instantiation

% start with a copy of the template
child = template;

% save global name for new instance
child.name = parent.name + "." + bind.localName;

fprintf("instantiating Base Component: %s\n",child.name);

% Apply variable mapping to all equations
% Traverse States & Transitions to reach all of them
numStates = length(child.States);
for st = 1:numStates
    St = child.States(st);
    
    % Apply mapping to Flow
    child.States(st).Flow = applyMappingToEquation(St.Flow,bind);
        
    % Apply Mapping to Invariant
    child.States(st).Invariant = applyMappingToCondition(St.Invariant,bind);
    
    numTrans = length(St.Trans);
    for tr = 1:numTrans
        Tr = St.Trans(tr);
        
        % Apply mapping to Reset Fct.
        child.States(st).Trans(tr).reset = applyMappingToEquation(Tr.reset,bind);
        
        % Apply mapping to Guard
        child.States(st).Trans(tr).guard = applyMappingToCondition(Tr.guard,bind);
    end
end


end