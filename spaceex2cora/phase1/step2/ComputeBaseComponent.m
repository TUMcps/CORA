function [bc_out,warnings] = ComputeBaseComponent( bc_in )
%INPUT:
%   bc_in (struct): component definition in SX format
%   (xml element converted to matlab struct)
%OUTPUT:
%   bc_out (struct): base component in structHA format
%   contains:listOfVar (variable definitions, cont. states/inputs/params)
%            States (list of discrete states (with flow & invariant))
%               (also includes transitions(with guard & reset))
%   CONVENTION: only strings leave this functions (no char arrays)
%   warnings (struct): document failed operations in warnings.messages
warnings = [];

% Assign meta information
bc_out.id = string(bc_in.Attributes.id);

% Collect Variables and Constants
[listOfVar, listOfLab] = CollectVariables(bc_in.param);
bc_out.listOfVar = listOfVar;
%bc_out.h_listOfLab = listOfLab; %unused, since CORA does not support this

% Compute States
% For each State we assign the meta information (name, id) and the
% transition array
num_states = length(bc_in.location);
for j = 1:num_states
    bc_out.States(j).id = string(bc_in.location{j}.Attributes.id);
    bc_out.States(j).name = string(bc_in.location{j}.Attributes.name);
    bc_out.States(j).Trans = [];
    
    % Assign the flow equation for this State
    if isfield(bc_in.location{j},'flow')
        % Retrieve equation text
        text = bc_in.location{j}.flow{1}.Text;
    else
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [vn,exprs,warn] = parseAssignment(text);
    % Store text & parsed equations in data structure
    bc_out.States(j).Flow.Text = string(text);
    bc_out.States(j).Flow.varNames = vn;
    bc_out.States(j).Flow.expressions = exprs;
    % store warnings
    warnings = [warnings,warn];
    
    
    % Assign the invariants for this State
    if isfield(bc_in.location{j},'invariant')
        %Retrieve equation text
        text = bc_in.location{j}.invariant{1}.Text;
    else
        % no equation => no conditions, global
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [ineqs,eqs,warn] = parseCondition(text);
    % Store text & parsed equations in data structure
    bc_out.States(j).Invariant.Text = string(text);
    bc_out.States(j).Invariant.inequalities = ineqs;
    bc_out.States(j).Invariant.equalities = eqs;
    % store warnings
    warnings = [warnings,warn];
end

% Compute number of Transitions, only for internal use.
if isfield(bc_in ,'transition')
    % If there is a fiel 'transition' then the length of the field gives
    % the number of transitions.
    h_numTrans = length(bc_in.transition);
else
    % If there is no field 'transition' then the number of transitions is
    % zero.
    h_numTrans = 0;
end

% For all transitions compute the source and the target of the transition
% to calculate the number of transitions for each state.
for k = 1:h_numTrans
    % Compute index of source & target state, for easier reference
    % resolution
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
    bc_out.States(h_source).Trans(h_trans).destination =  h_target;
    
    
    % Assign guards if the field 'guard' is given.
    if isfield(bc_in.transition{k},'guard')
        %Retrieve equation text
        text = bc_in.transition{k}.guard{1}.Text;
    else
        % no equation => no conditions, global
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [ineqs,eqs,warn] = parseCondition(text);
    % Store text & parsed equations in data structure
    bc_out.States(h_source).Trans(h_trans).guard.Text = string(text);
    bc_out.States(h_source).Trans(h_trans).guard.inequalities = ineqs;
    bc_out.States(h_source).Trans(h_trans).guard.equalities = eqs;
    % store warnings
    warnings = [warnings,warn];
    
    % Assign jump function if the field 'assignment' is given.
    if isfield(bc_in.transition{k},'assignment')
        % Retrieve equation text
        text = bc_in.transition{k}.assignment{1}.Text;
    else
        % no equation => reset is identity function
        text = "";
    end
    % Parse text into a partially symbolic form, to simplify variable
    % substitution
    [vn,exprs,warn] = parseAssignment(text);
    % Store text & parsed equations in data structure
    bc_out.States(h_source).Trans(h_trans).reset.Text = string(text);
    bc_out.States(h_source).Trans(h_trans).reset.varNames = vn;
    bc_out.States(h_source).Trans(h_trans).reset.expressions = exprs;
    % store warnings
    warnings = [warnings,warn];
    
end

end


