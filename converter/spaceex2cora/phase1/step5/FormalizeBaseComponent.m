function BC = FormalizeBaseComponent(BC,convtype)
% FormalizeBaseComponent - splits variables into inputs, states, constants;
%    applies equation-parsing scripts to fields given as string equations
%
% Description:
%    eq2linSys to:     BC.States{*}.Flow,
%                      BC.States{*}.Trans{*}.reset
%    eq2polytope to:   BC.States{*}.Invariant
%                      BC.States{*}.Trans{*}.guard
%
% Syntax:
%    BC = FormalizeBaseComponent(comp,convtype)
%
% Inputs:
%    BC (struct) - base component to be formalized
%    convtype - true if parallel hybrid automaton to be generated, false
%               if flat hybrid automaton to be generated
%
% Outputs:
%    BC - formalized based component containing the flow equation,
%         invariant, guard set, and reset function as matrices
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   22-June-2022 (MW, empty sets for invariants, guards)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Unfortunately, CORA does not support 0-input systems yet.
% If the system does not have inputs, we add a dummy input without effect.
if isempty(BC.inputs)
    % generate a dummy input "uDummy" (with as many underscores at the end
    % as necessary until it is unique; in principle uDummy could already be
    % used in the model)
    maxlength = 0;
    for i = 1:length(BC.listOfVar)
        if regexp(BC.listOfVar(i).name,'^uDummy(_*)$')
            maxlength = max(maxlength,strlength(BC.listOfVar(i).name));
        end
    end
    % build dummy name as char vector, add underscores until it is unique
    dummyName = 'uDummy';
    dummyName = [dummyName repmat('_',1,maxlength+1-length(dummyName))];
    
    % convert name to string and add dummy input variable
    BC.inputs(1).name = string(dummyName);
end

% iterate over locations of base component
for i = 1:length(BC.States)
    % read out i-th location for easier access
    State = BC.States(i);
    
    % derive linear representation of flow equations
    if convtype
        [isLin,A,B,c,eqs] = eq2linSys(State.Flow.varNames,...
            State.Flow.expressions,BC.states,BC.inputs);
        [isLin_out,C,D,k,eqs_out] = outputMatricesComp2Comp(BC.states,BC.inputs,...
            [BC.outputsLocal;BC.outputsGlobal],State.Invariant);
    else
        % what do outputsLocal map to in state vector
        if ~isempty(BC.outputsLocal)
            map = string(State.Invariant.exprRight(...
                strcmp(string(State.Invariant.exprLeft),BC.outputsLocal.name)));
        else
            map = "";
        end
        % derive linear representation of flow equations
        [isLin,A,B,c,eqs] = eq2linSysFlat(State.Flow.varNames,...
            State.Flow.expressions,BC.states,BC.inputs,BC.outputsLocal,map);
        % output matrices
        [isLin_out,C,D,k,eqs_out] = outputMatrices(State.Invariant.equalities,...
            BC.states,BC.inputs,BC.outputsLocal,map,BC.outputsGlobal);
        % write output equation to invariant for later
        BC.States(i).Invariant.Text_output = eqs_out;
    end

    % assign matrices if system is linear
    if isLin && isLin_out
        BC.States(i).Flow.A = A;
        BC.States(i).Flow.B = B;
        BC.States(i).Flow.c = c;
        % output matrices relate local outputs to inputs of other
        % components (only parallel hybrid automata)
        BC.States(i).Flow.C = C;
        BC.States(i).Flow.D = D;
        BC.States(i).Flow.k = k;
    end
    % save differential equations in format: dx(i,1) = ...
    BC.States(i).Flow.FormalEqs = eqs;
    % save output equations in format: dx(i,1) = ... (we can write x
    % instead of y since the computation is outsourced to a function)
    BC.States(i).Flow.FormalEqs_out = eqs_out;
        
    % derive set for invariant
    if State.Invariant.Text == ""
        % empty invariant: short version to define fullspace as invariant
        set = fullspace(length(State.Flow.expressions));
    else
        % non-empty invariant
        set = eq2set(State.Invariant.inequalities,State.Invariant.equalities,...
                 BC.states,[BC.outputsLocal BC.outputsGlobal]);
    end
    BC.States(i).Invariant.set = set;
    
    % iterate over outgoing transitions of current State
    if isfield(State,'Trans')
        nrTrans = length(State.Trans);
    else
        nrTrans = 0;
    end

    for j = 1:nrTrans
        Tran = State.Trans(j);
        
        %derive polytope for guard
        if strcmp(Tran.guard.Text,"")
            % no guard set given -> instant transition (fullspace)
            set = fullspace(length(State.Flow.expressions));
        elseif convtype
            % guard set given, pHA
            set = eq2set(Tran.guard.inequalities,Tran.guard.equalities, ...
                         BC.states,BC.outputsGlobal);
        else
            % guard set given, flat HA
            set = eq2set(Tran.guard.inequalities,Tran.guard.equalities,...
                         BC.states,[BC.outputsLocal BC.outputsGlobal]);
        end
        BC.States(i).Trans(j).guard.set = set;
        
        % derive linear representation of assignment
        [isLin,A,B,c,eqs,containsInput] = eq2linSys(Tran.reset.varNames,...
            Tran.reset.expressions,BC.states,BC.inputs,'assignment');
        if isLin
            BC.States(i).Trans(j).reset.A = A;
            BC.States(i).Trans(j).reset.c = c;
            if containsInput
                % reset is input-dependent
                BC.States(i).Trans(j).reset.B = B;
            end
        else
            if strlength(eqs) ~= 0
                BC.States(i).Trans(j).reset.FormalEqs = eqs;
            else
                warning("Transition from ""%s"" to ""%s"" does not have a reset function. " + ...
                    "Using identity...",State.name,BC.States(Tran.destination).name);
                BC.States(i).Trans(j).reset.A = diag(ones(size(BC.states)));
                BC.States(i).Trans(j).reset.c = zeros(size(BC.states));
            end
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
