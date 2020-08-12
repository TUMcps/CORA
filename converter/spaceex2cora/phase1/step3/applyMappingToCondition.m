function [condition] = applyMappingToCondition(condition,bind)
% Applies the variable mapping of a Bind struct to a Invariant or Guard struct
% INPUTS
%   condition (struct): Invariant or Guard struct in BaseComponent format
%                       (needed fields: inequalities,equalities)
%   bind (struct): Bind struct in NetworkComponent format
%                       (needed fields: keys,values,renames)
% OUTPUTS
%   condition with variable mappings applied

% substitute variables for mapped values in Inequations
condition.inequalities = applySymMapping(condition.inequalities,bind.keys,bind.values);
% substitute variables for mapped values in Equations
condition.equalities = applySymMapping(condition.equalities,bind.keys,bind.values);
end

