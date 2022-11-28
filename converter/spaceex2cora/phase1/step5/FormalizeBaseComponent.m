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

% Author:       ???
% Written:      ???
% Last update:  22-June-2022 (MW, empty sets for invariants, guards)
% Last revision:---

%------------- BEGIN CODE --------------

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

% iterate over discrete States
for i = 1:length(BC.States)
    State = BC.States(i);
    
    % derive linear representation of flow equations
    if convtype
        [isLin,A,B,c,eqs] = eq2linSys(State.Flow.varNames,...
            State.Flow.expressions,BC.states,BC.inputs);
        [C,D,k] = outputMatricesComp2Comp(BC.states,BC.inputs,BC.outputsLocal);
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
    end

    % assign matrices for if system is linear
    if isLin
        BC.States(i).Flow.A = A;
        BC.States(i).Flow.B = B;
        BC.States(i).Flow.c = c;
        if convtype
            % only parallel: assign output matrices (to be used as inputs
            % in other components)
            BC.States(i).Flow.C = C;
            BC.States(i).Flow.D = D;
            BC.States(i).Flow.k = k;
        end
    end
    % save equations in format: dx(i,1) = ...
    BC.States(i).Flow.FormalEqs = eqs;
        
    % derive polytope for invariant
    if State.Invariant.Text == ""
        % empty invariant: short version to define full space as invariant
        set = [];
    else
        % non-empty invariant
        set = eq2set(State.Invariant.inequalities,State.Invariant.equalities,...
                 BC.states,[BC.outputsLocal BC.outputsGlobal]);
    end
    BC.States(i).Invariant.set = set;
    
    % assign linear output matrices: only flat HA
    % note: this might be overhauled in the future as it is not clear
    % what is the relation between the invariant and the outputs
    if ~convtype
        [C,D,k] = outputMatrices(State.Invariant.equalities,...
            BC.states,BC.inputs,BC.outputsGlobal);
        BC.States(i).Flow.C = C;
        BC.States(i).Flow.D = D;
        BC.States(i).Flow.k = k;
    end
    
    % iterate over outgoing transitions of current State
    if isfield(State,'Trans')
        numTrans = length(State.Trans);
    else
        numTrans = 0;
    end
    for j = 1:numTrans
        Tran = State.Trans(j);
        
        %derive polytope for guard
        if strcmp(Tran.guard.Text,"")
            % no guard set given -> immediate transition (guard = [])
            set = [];
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

%------------- END OF CODE --------------
