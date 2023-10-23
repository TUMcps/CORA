function [bc_out,warnings] = ComputeBaseComponent(bc_in)
% ComputeBaseComponent - For a base component of the SpaceEx model, extract
%    properties of the component (variable names, expressions for
%    flows/guards/resets/etc.) from the struct and prepare them for later
%    conversion to CORA format;
%    note: only strings leave this function (no char arrays)!
%
% Syntax:
%    [bc_out,warnings] = ComputeBaseComponent(bc_in)
%
% Inputs:
%    bc_in (struct) - component definition in SX format
%                     (part of SpaceEx xml-file converted to struct)
%
% Outputs:
%    bc_out (struct) - base component in structHA format, with fields
%            .listOfVar (variable definitions, cont. states/inputs/params)
%            .States (list of discrete states (with flow & invariant),
%                     and transitions (with guard & reset))
%    warnings (struct) - document failed operations in warnings.messages
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       ---
% Last update:   12-January-2023 (MW, add read-out of output equations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init warnings
warnings = [];

% Assign meta-data
bc_out.id = string(bc_in.Attributes.id);

% Collect Variables and Constants
[listOfVars, listOfLabels] = CollectVariables(bc_in.param);
bc_out.listOfVar = listOfVars;
bc_out.listOfLabels = listOfLabels;

% go over each discrete state aka location: assign
% - the meta-data (name, id)
% - the transition array (default: [])
% - the flow equation
% - the invariant
num_states = length(bc_in.location);
for j = 1:num_states

    % meta-data and initialization of outgoing transitions
    bc_out.States(j).id = string(bc_in.location{j}.Attributes.id);
    bc_out.States(j).name = string(bc_in.location{j}.Attributes.name);
    bc_out.States(j).Trans = [];
    
    % flow equation
    if isfield(bc_in.location{j},'flow')
        % Retrieve equation text
        text = bc_in.location{j}.flow{1}.Text;
    else
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [vn,exprs,warn] = parseAssignment(text);
    % store text (as written in SpaceEx editor)
    bc_out.States(j).Flow.Text = string(text);
    % store variable names (left-hand side of flow expressions)
    bc_out.States(j).Flow.varNames = vn;
    % store derivatives (right-hand side of flow expressions)
    bc_out.States(j).Flow.expressions = exprs;
    % store warnings
    warnings = [warnings,warn];
    

    % invariant
    if isfield(bc_in.location{j},'invariant')
        %Retrieve equation text
        text = bc_in.location{j}.invariant{1}.Text;
    else
        % no equation => no conditions, global
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [ineqs,eqs,expLeft,expRight,warn] = parseCondition(text);
    % store text (as written in SpaceEx editor)
    bc_out.States(j).Invariant.Text = string(text);
    % store text of output equations
    bc_out.States(j).Invariant.Text_output = ...
        aux_readOutputEqFromInvariant(string(text));
    % store symbolic inequalities
    bc_out.States(j).Invariant.inequalities = ineqs;
    % store symbolic equalities
    bc_out.States(j).Invariant.equalities = eqs;
    % symbolic variables occurring on left- and right-hand side
    % (this will be relevant for the conversion to flat hybrid automata)
    bc_out.States(j).Invariant.exprLeft = expLeft;
    bc_out.States(j).Invariant.exprRight = expRight;
    % store warnings
    warnings = [warnings,warn];
end

% number of transitions
h_numTrans = 0;
if isfield(bc_in,'transition')
    % length of the field gives the number of transitions
    h_numTrans = length(bc_in.transition);
end

% loop over all transitions, read
% - source and target location of the transition
% - guard set
% - reset function
% - number of transitions for each location
for k = 1:h_numTrans

    % source and target location
    source_id = bc_in.transition{k}.Attributes.source;
    target_id = bc_in.transition{k}.Attributes.target;
    h_source = 0;
    h_target = 0;
    for j = 1:num_states
        if strcmp(bc_out.States(j).id,source_id)
            h_source = j;
        end
        if strcmp(bc_out.States(j).id,target_id)
            h_target = j;
        end
    end
    % Increment number of transitions of State h_source, save in h_trans
    h_trans = length(bc_out.States(h_source).Trans) + 1;
    
    % Assign destination
    bc_out.States(h_source).Trans(h_trans).destination = h_target;
    
    
    % guard set
    if isfield(bc_in.transition{k},'guard')
        % retrieve equation text
        text = bc_in.transition{k}.guard{1}.Text;
    else
        % no equation -> guard set is entire space (transition is taken
        % instantly upon entering the location unless synchronization
        % label is given)
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [ineqs,eqs,~,~,warn] = parseCondition(text);
    % store text (as written in SpaceEx editor)
    bc_out.States(h_source).Trans(h_trans).guard.Text = string(text);
    % store symbolic inequalities of guard set
    bc_out.States(h_source).Trans(h_trans).guard.inequalities = ineqs;
    % store symbolic equalities of guard set
    bc_out.States(h_source).Trans(h_trans).guard.equalities = eqs;
    % store warnings
    warnings = [warnings,warn];
    

    % reset function
    if isfield(bc_in.transition{k},'assignment')
        % retrieve equation text
        text = bc_in.transition{k}.assignment{1}.Text;
    else
        % no equation -> reset is identity function (to be processed later)
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [vn,exprs,warn] = parseAssignment(text);
    % store text (as written in SpaceEx editor)
    bc_out.States(h_source).Trans(h_trans).reset.Text = string(text);
    % store symbolic variable names (left-hand side of reset functions)
    bc_out.States(h_source).Trans(h_trans).reset.varNames = vn;
    % store symbolic expressions (right-hand side of reset functions)
    bc_out.States(h_source).Trans(h_trans).reset.expressions = exprs;
    % store warnings
    warnings = [warnings,warn];
    
    
    % store synchronization labels
    if isfield(bc_in.transition{k},'label')
        bc_out.States(h_source).Trans(h_trans).label = ...
            bc_in.transition{k}.label{1}.Text;
    else
        % no label => transition is not synchronized with any other
        bc_out.States(h_source).Trans(h_trans).label = "";
    end
    
end

end


% Auxiliary functions -----------------------------------------------------

function outputEqStr = aux_readOutputEqFromInvariant(str)

% check whether there are any output equations at all
if ~contains(str,"==")
    outputEqStr = ""; return
end

% only one equation
if ~contains(str,"&&")
    outputEqStr = str; return
end

% we have some (actual) invariants with comparison operators other than
% "==" and output equations using the comparison operator "=="

% indices indicating output equations
idxEq = strfind(str,"==");
% indices where equations are separated
idxAnd = strfind(str,"&&");

% container for output equations
outputEqs = {};
% go over equations one by one
nrEqs = length(idxAnd)+1;
for i=1:nrEqs
    if i == 1
        eq = extractBefore(str,idxAnd(1));
        % first equation is an output equation
        if contains(eq,"==")
            outputEqs = [outputEqs; {eq}];
        end
    elseif i < nrEqs
        eq = extractAfter(str,idxAnd(i-1)+1);
        shift = strlength(str) - strlength(eq);
        eq = extractBefore(eq,idxAnd(i)-shift);
        if contains(eq,"==")
            outputEqs = [outputEqs; {eq}];
        end
    else % last equation
        eq = extractAfter(str,idxAnd(end)+1);
        if contains(eq,"==")
            outputEqs = [outputEqs; {eq}];
        end
    end
end

% put equations into nicely readable form (note: due to the removal of
% invariant inequalities, the original format has already been changed)
for i=1:length(outputEqs)
    % truncate whitespaces at the beginning and at the end
    while true
        % logical array of whitespaces
        whitespaces = isstrprop(outputEqs{i},'wspace');
        if ~whitespaces(1)
            break
        else
            outputEqs{i} = extractAfter(outputEqs{i},1);
        end
    end
    while true
        % logical array of whitespaces
        whitespaces = isstrprop(outputEqs{i},'wspace');
        if ~whitespaces(end)
            break
        else
            outputEqs{i} = extractBefore(outputEqs{i},strlength(outputEqs{i}));
        end
    end
    % convert to char to enable call of strjoin with a cell-array
    outputEqs{i} = char(outputEqs{i});
end
outputEqStr = string(strjoin(outputEqs," && "));

end

% ------------------------------ END OF CODE ------------------------------
