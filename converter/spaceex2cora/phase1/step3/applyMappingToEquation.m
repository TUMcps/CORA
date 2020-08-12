function [equation] = applyMappingToEquation(equation,bind)
% Applies the variable mapping of a Bind struct to a Flow or Reset struct
% INPUTS
%   equation (struct): Flow or Reset struct in BaseComponent format
%                       (needed fields: varnames,exprs)
%   bind (struct): Bind struct in NetworkComponent format
%                       (needed fields: keys,values,renames)
% OUTPUTS
%   equation with variable mappings applied

% rename left-side variables
equation.varNames = applyRenames(equation.varNames,bind.keys,bind.renames);
% substitute variables for mapped values in right-side expressions
equation.expressions = applySymMapping(equation.expressions,bind.keys,bind.values);
end

