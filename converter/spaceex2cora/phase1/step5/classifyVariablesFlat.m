function BC = classifyVariablesFlat(BC)
% classifyVariablesFlat - variables are classified as states, inputs, or
%    outputs, depending on where they appear in the given flow equations
%    note: only for flat hybrid automata
%
% Syntax:
%    BC = classifyVariablesFlat(BC)
%
% Inputs:
%    BC - merged base component resulting from automaton product
%
% Outputs:
%    BC (struct) - base component extended by list of variables, with
%       .states - variables appearing on the left side
%       .inputs - variables appearing only on the right side
%       .outputsLocal - output variables which are inputs to other loc
%                       (only occurring in originally parallel structures)
%       .outputsGlobal - output variables of the entire flat automaton
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: classifyVariablesParallel
%
% References: 
%   -

% Authors:       Mark Wetzlinger
% Written:       25-June-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% names of states (left-hand side in flow equations)
allStateNames = [];
% right-hand side of flow equations
allFlowExprs = [];

for i = 1:length(BC.States)
    % reminder: assignment expressions are stored in column vectors
    allStateNames = [allStateNames; BC.States(i).Flow.varNames];
    allFlowExprs = [allFlowExprs; BC.States(i).Flow.expressions];
    
    for j = 1:length(BC.States(i).Trans)
        allStateNames = [allStateNames; BC.States(i).Trans(j).reset.varNames];
        % If inputs/constants are ever changed within a reset (not sure if
        % possible), this would lead to an error; in this case, make a new
        % variable "resetExpressions" and integrate it into classifyVariables.m
        allFlowExprs = [allFlowExprs; BC.States(i).Trans(j).reset.expressions];
    end
end
% remove duplicates
allStateNames = unique(allStateNames);
allFlowExprs = unique(allFlowExprs);

if isempty(allStateNames)
    warning("No flow equations given; assumption: " + ...
        "all states occur in location invariants!");
    for i=1:length(BC.States)
        names = symvar([BC.States(i).Invariant.equalities;...
            BC.States(i).Invariant.inequalities]);
        allStateNames = [allStateNames;string(names)'];
    end
end

% difference in conversion to flatHA / parallelHA
invExprsLeft = [];
invExprsRight = [];
for i=1:length(BC.States)
    invExprsLeft = [invExprsLeft; BC.States(i).Invariant.exprLeft];
    invExprsRight = [invExprsRight; BC.States(i).Invariant.exprRight];
end

% number of variables
num_vars = length(BC.listOfVar);

% extract symbolic variables from right equation sides
% symvar returns the union set, if multiple symbolics are passed
flow_expr_vars = symvar(allFlowExprs);
% convert to variable names
flow_expr_varnames = string(flow_expr_vars);

% check whether variables appear in left flow sides (varnames)
flowleftIdx = false(1,num_vars);
% check whether variables appear in right flow sides (flowexprs)
flowrightIdx = false(1,num_vars);

for i=1:num_vars
    flowleftIdx(i) = any(allStateNames == BC.listOfVar(i).name);
    flowrightIdx(i) = any(flow_expr_varnames == BC.listOfVar(i).name);
end

% extension for flat conversion
inv_expr_left_vars = symvar(invExprsLeft);
inv_expr_right_vars = symvar(invExprsRight);
inv_expr_left_varnames = string(inv_expr_left_vars);
inv_expr_right_varnames = string(inv_expr_right_vars);
% check whether variables appear in left inv sides (invexprs)
invleftIdx = false(1,num_vars);
% check whether variables appear in right inv sides (invexprs)
invrightIdx = false(1,num_vars);
for i = 1:num_vars
    invleftIdx(i) = any(inv_expr_left_varnames == BC.listOfVar(i).name);
    invrightIdx(i) = any(inv_expr_right_varnames == BC.listOfVar(i).name);
end

% chop up listOfVar using logical indexing

% 1. states appear on the left-hand side of the flow equation
BC.states = BC.listOfVar(flowleftIdx);

% 2. inputs appear on the right-hand side of the flow equation, but not on
% the left-hand side, additionally, not on the left-hand side of the
% invariant
BC.inputs = BC.listOfVar(~flowleftIdx & flowrightIdx & ~invleftIdx);

% 3. local outputs: outputs in what was previously a parallel base
% component, but is now part of the merged component
BC.outputsLocal = BC.listOfVar(flowrightIdx & invleftIdx);

% 4. global outputs: outputs of the entire automaton
BC.outputsGlobal = BC.listOfVar(~flowrightIdx & invleftIdx);

% ------------------------------ END OF CODE ------------------------------
