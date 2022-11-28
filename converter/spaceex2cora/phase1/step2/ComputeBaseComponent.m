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

% Author:       ???
% Written:      ???
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init warnings
warnings = [];

% Assign meta-data
bc_out.id = string(bc_in.Attributes.id);

% Collect Variables and Constants
[listOfVars, listOfLabels] = CollectVariables(bc_in.param);
bc_out.listOfVar = listOfVars;
%bc_out.listOfLabels = listOfLabels; % unused

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
        % immediately upon entering the location unless synchronization
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
    
    
    % synchronization labels
    if isfield(bc_in.transition{k},'label')
        syncLabel = bc_in.transition{k}.label{1}.Text;
    else
        % no label => transition is not synchronized
        syncLabel = "";
    end
    % Store label in data structure
    bc_out.States(h_source).Trans(h_trans).label = string(syncLabel);
    
end

%------------- END OF CODE --------------
