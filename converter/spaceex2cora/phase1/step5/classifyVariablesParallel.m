function BClist = classifyVariablesParallel(BClist)
% classifyVariablesParallel - variables are classified as states, inputs,
%    or outputs, depending on where they appear in the given flow equations
%
% Syntax:
%    BClist = classifyVariablesParallel(BClist)
%
% Inputs:
%    BClist - list of base components
%
% Outputs:
%    BClist - list of base components extended by 
%       .states - variables appearing on the left side
%       .inputs - variables appearing only on the right side
%       .outputsLocal - output variables which are inputs to other loc
%       .outputsGlobal - output variables of the entire parallel automaton
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: classifyVariablesFlat
%
% References: 
%   -

% Author:        Mark Wetzlinger
% Written:       25-June-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

% number of base components
numBC = length(BClist);

% initialize variables for state names (left-hand side of flow equation)
% and flow expressions (right-hand side of flow equation)
allStateNames = cell(numBC,1);
allFlowExprs = cell(numBC,1);

% loop over all base components
for i=1:length(BClist)
    
    for j = 1:length(BClist(i).States)
        % reminder: assignment expressions are stored in column vectors
        allStateNames{i} = [allStateNames{i}; BClist(i).States(j).Flow.varNames];
        allFlowExprs{i} = [allFlowExprs{i}; BClist(i).States(j).Flow.expressions];
        
        for k = 1:length(BClist(i).States(j).Trans)
            % add left-hand side of outgoing reset equations 
            % -> shouldn't this be the incoming reset equations?
            allStateNames{i} = [allStateNames{i}; BClist(i).States(j).Trans(k).reset.varNames];
            % If inputs/constants are ever changed within a reset (not sure if
            % possible), this would lead to an error; in this case, make a new
            % variable "resetExpressions" and integrate it into classifyVariables.m
            allFlowExprs{i} = [allFlowExprs{i}; BClist(i).States(j).Trans(k).reset.expressions];
        end
    end
    
    % fallback solution in case no location has a flow equation
    % (users should, however, not rely too much on this feature; it is
    % better to properly define all flow equations)
    if isempty(allStateNames{i})
        warning("No flow equations given; assumption: all states occur in location invariants!");
        for j = 1:length(BClist(i).States)
            names = symvar([BClist(i).States(j).Invariant.equalities;...
                BClist(i).States(j).Invariant.inequalities]);
            allStateNames{i} = [allStateNames{i}; string(names)'];
        end
    end

    % number of variables
    num_vars = length(BClist(i).listOfVar);
    
    % extract symbolic variables from right equation sides
    % symvar returns the union set, if multiple symbolics are passed
    flow_expr_varnames = string(symvar(allFlowExprs{i}));
    
    % check whether variables appear in left flow sides (varnames)
    flowleftIdx = false(1,num_vars);
    % check whether variables appear in right flow sides (flowexprs)
    flowrightIdx = false(1,num_vars);
    
    for j=1:num_vars
        flowleftIdx(j) = any(allStateNames{i} == BClist(i).listOfVar(j).name);
        flowrightIdx(j) = any(flow_expr_varnames == BClist(i).listOfVar(j).name);
    end
    
    % states appear on the left-hand side of the flow equation
    BClist(i).states = BClist(i).listOfVar(flowleftIdx);
    
    % inputs appear on the right-hand side of the flow equation, but not on
    % the left-hand side
    BClist(i).inputs = BClist(i).listOfVar(~flowleftIdx & flowrightIdx);

    % init local outputs = inputs in other components
    % (determined in next loop)
    BClist(i).outputsLocal = BClist(i).listOfVar(false(num_vars,1));

    % no global outputs, since SpaceEx does not support global outputs
    BClist(i).outputsGlobal = [];
end

% outputs require knowledge about states, inputs in all base components,
% thus we can only now determine the outputs for each component
for i=1:numBC
    % loop over all other components and see whether they use a state
    % variable of this component as an input
    allOutputNames = [];
    varnames = [BClist(i).listOfVar.name];
    for j=1:numBC
        if i == j || isempty(BClist(j).inputs)
            % skip self-reading and empty input list of other component
            continue;
        end
        % which members of the list of variables are inputs in component j?
        outputIdx = ismember(varnames,[BClist(j).inputs.name]);
        % catch repetitions: do not use unique (after the loop) in order to
        % preserve the order of the variables (should give nicer output
        % matrices than when alphabetical order is used)
        for k=1:length(outputIdx)
            if outputIdx(k) && ~isempty(allOutputNames)
                if any(allOutputNames == BClist(j).inputs(k).name)
                    outputIdx(k) = false;
                end
            end
        end
        allOutputNames = [allOutputNames; varnames(outputIdx)'];
    end
    for j=1:length(allOutputNames)
        BClist(i).outputsLocal(j).name = allOutputNames(j);
    end

end

%------------- END OF CODE -------------
