function [states,inputs,outputs] = classifyVariables(listOfVar,varnames,exprs)
%	INTPUTS:
%   listOfVar (struct): variables of a component
%           (output of CollectVaribles)
%   varnames/exprs (string/symbolic): parsed flow equation
%           (output of parseAssignment)
%           to classify based on mutiple flows, just concat their
%           varnames & exprs arrays (order of entries is irrelevant)
%   OUTPUTS:
%   returns a partition of the input listOfVar
%   variables are classified as states,inputs or outputs, depending on
%   where they appear in the given flow equations
%   states: variables appearing on the left side
%   inputs: variables appearing only on the right side
%   outputs: variables appearing nowhere
num_vars = length(listOfVar);

% extract symbolic variables from right equation sides
% symvar returns the union set, if multiple symbolics are passed
expr_vars = symvar(exprs);
% convert to variable names
expr_varnames = string(expr_vars);

% check, wether variables appear in left sides (varnames)
leftIdx = false(1,num_vars);
for i = 1:num_vars
    leftIdx(i) = any(varnames == listOfVar(i).name);
end

% check, wether variables appear in right sides (exprs)
rightIdx = false(1,num_vars);
for i = 1:num_vars
    rightIdx(i) = any(expr_varnames == listOfVar(i).name);
end

% chop up listOfVar using logical indexing
states = listOfVar(leftIdx);
inputs = listOfVar(and(~leftIdx,rightIdx));
outputs = listOfVar(and(~leftIdx,~rightIdx));

end

