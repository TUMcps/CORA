function [states,inputs,outputsLocal,outputsGlobal] = ...
    classifyVariables(listOfVar,varnames,flowexprs,invexprsleft,invexprsright)
% classifyVariables - variables are classified as states,
%   inputs or outputs, depending on where they appear in the given flow equations
%
% Syntax:  
%    [states,inputs,outputsLocal,outputsGlobal] = ...
%      classifyVariables(listOfVar,varnames,flowexprs,invexprsleft,invexprsright)
%
% Inputs:
%   listOfVar     - variables of a component (output of CollectVariables)
%   varnames      - parsed flow equation (output of parseAssignment)
%           to classify based on multiple flows, just concat their
%           varnames & exprs arrays (order of entries is irrelevant)
%   flowexprs     - 
%   invexprsleft  - 
%   invexprsright - 
%
% Outputs:
%   states        - variables appearing on the left side
%   inputs        - variables appearing only on the right side
%   outputsLocal  - output variables which are inputs to other loc
%   outputsGlobal - global output variables
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: 
%   -

% Author:        Mark Wetzlinger
% Written:       ??-??-2019
% Last update:   26-February-2020
% Last revision: ---

%------------- BEGIN CODE --------------

if nargin == 3 || nargin == 5
    % ok
else
    error("Function classifyVariables has an invalid number of input parameters");
end

num_vars = length(listOfVar);

% extract symbolic variables from right equation sides
% symvar returns the union set, if multiple symbolics are passed
flow_expr_vars = symvar(flowexprs);
% convert to variable names
flow_expr_varnames = string(flow_expr_vars);

% check whether variables appear in left flow sides (varnames)
flowleftIdx = false(1,num_vars);
% check whether variables appear in right flow sides (flowexprs)
flowrightIdx = false(1,num_vars);

for i = 1:num_vars
    flowleftIdx(i) = any(varnames == listOfVar(i).name);
    flowrightIdx(i) = any(flow_expr_varnames == listOfVar(i).name);
end

% extension for flat conversion
if nargin == 5
    inv_expr_left_vars = symvar(invexprsleft);
    inv_expr_right_vars = symvar(invexprsright);
    inv_expr_left_varnames = string(inv_expr_left_vars);
    inv_expr_right_varnames = string(inv_expr_right_vars);
    % check whether variables appear in left inv sides (invexprs)
    invleftIdx = false(1,num_vars);
    % check whether variables appear in right inv sides (invexprs)
    invrightIdx = false(1,num_vars);
    
    for i = 1:num_vars
        invleftIdx(i) = any(inv_expr_left_varnames == listOfVar(i).name);
        invrightIdx(i) = any(inv_expr_right_varnames == listOfVar(i).name);
    end
end

% chop up listOfVar using logical indexing
states = listOfVar(flowleftIdx);

if nargin == 3
    inputs = listOfVar(and(~flowleftIdx,flowrightIdx));
    % outputsGlobal since component doesn't use its own output
    outputsGlobal = listOfVar(and(~flowleftIdx,~flowrightIdx));
    outputsLocal = [];
else
    inputs = listOfVar(and(and(~flowleftIdx,flowrightIdx),~invleftIdx));
    % flat merged automaton has local and global outputs
    outputsLocal = listOfVar(and(flowrightIdx,invleftIdx));
    outputsGlobal = listOfVar(and(~flowrightIdx,invleftIdx));
end

end

