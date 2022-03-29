function BC = FormalizeBaseComponent(comp,convtype)
% splits variables into inputs,states & constants
% applies equation-parsing scripts to fields given as string equations
% eq2linSys to:     BC.States{*}.Flow,
%                   BC.States{*}.Trans{*}.reset
% eq2polytope to:   BC.States{*}.Invariant
%                   BC.States{*}.Trans{*}.guard
% convtype: parallelHA - 1; flatHA - 0

BC = comp;

% split variables into states, inputs & others
% this requires gathering the flow equations of all states
% as well as the reset expressions of all transitions
allStateNames = [];
allFlowExprs = [];

for i = 1:length(BC.States)
    % reminder: assignment expressions are stored in column vectors
    allStateNames = [allStateNames; BC.States(i).Flow.varNames];
    
    allFlowExprs = [allFlowExprs; BC.States(i).Flow.expressions];
    
    for j = 1:length(BC.States(i).Trans)
        allStateNames = [allStateNames; BC.States(i).Trans(j).reset.varNames];
        % If inputs/constants are ever changed within a reset (not sure if possible)
        % this will lead to an error, in this case, make new variable
        % "resetExpressions" and integrate it into classifyVariables.m
        allFlowExprs = [allFlowExprs; BC.States(i).Trans(j).reset.expressions];
    end
end

if isempty(allStateNames)
    warning("No flow equations given, assuming all states to be occuring within location invariants!");
    for i = 1:length(BC.States)
        names = symvar([BC.States(i).Invariant.equalities;BC.States(i).Invariant.inequalities]);
        allStateNames = [allStateNames;string(names)'];
    end
end

% difference in conversion to flatHA / parallelHA
if ~convtype
    % flatHA:
    InvExprsLeft = [];
    InvExprsRight = [];
    for i = 1:length(BC.States)
        InvExprsLeft = [InvExprsLeft; BC.States(i).Invariant.exprLeft];
        InvExprsRight = [InvExprsRight; BC.States(i).Invariant.exprRight];
    end
    [BC.states,BC.inputs,BC.outputsLocal,BC.outputsGlobal] = ...
        classifyVariables(BC.listOfVar,allStateNames,allFlowExprs,InvExprsLeft,InvExprsRight);
else
    % parallelHA / (non-)linear:
    [BC.states,BC.inputs,BC.outputsLocal,BC.outputsGlobal] = ...
        classifyVariables(BC.listOfVar,allStateNames,allFlowExprs);
end

% variables that are assigned a derivative are states
% variables that influence a derivative are inputs/parameters
% other variables are currently being ignored

% Unfortunately, CORA does not support 0-input systems yet.
% If the system is input-less, add a dummy input without effect.
if isempty(BC.inputs)
    % generate a dummy input "uDummy"
    % (make sure dummy name is not used as a variable)
    maxlength = 0;
    for i = 1:length(BC.listOfVar)
        if regexp(BC.listOfVar(i).name,'^uDummy(_*)$')
            maxlength = max(maxlength,strlength(BC.listOfVar(i).name));
        end
    end
    % build dummy name as char vector, add understrikes until it is unique
    dummyName = 'uDummy';
    dummyName = [dummyName repmat('_',1,maxlength+1-length(dummyName))];
    % this may be sliiiiight overkill, but it is 100% safe
    
    % convert name to string & add dummy input variable
    BC.inputs(1).name = string(dummyName);
end

% iterate over discrete States
for i = 1:length(BC.States)
    State = BC.States(i);
    
    % derive linear representation of flow equations
    if convtype
        [isLin,A,B,c,eqs] = eq2linSys(State.Flow.varNames,State.Flow.expressions,BC.states,BC.inputs);
    else
        % what do outputsLocal map to in state vector
        if ~isempty(BC.outputsLocal)
            map = string(State.Invariant.exprRight(strcmp(string(State.Invariant.exprLeft),BC.outputsLocal.name)));
        else
            map = "";
        end
        % derive linear representation of flow equations
        [isLin,A,B,c,eqs] = eq2linSysFlat(State.Flow.varNames,State.Flow.expressions,...
            BC.states,BC.inputs,BC.outputsLocal,map);
    end
    if(isLin)
        BC.States(i).Flow.A = A;
        BC.States(i).Flow.B = B;
        BC.States(i).Flow.c = c;
    end
    
    BC.States(i).Flow.FormalEqs = eqs;
        
    % derive polytope for invariant
    %[A,b,Ae,be] = eq2polytope(State.Invariant.inequalities,State.Invariant.equalities,BC.states,BC.outputs);
    set = eq2set(State.Invariant.inequalities,State.Invariant.equalities,...
                 BC.states,[BC.outputsLocal BC.outputsGlobal]);
    BC.States(i).Invariant.set = set;
    
    [C,D,k] = outputMatrices(State.Invariant.equalities,BC.states,BC.inputs,BC.outputsGlobal);
    BC.States(i).Invariant.C = C;
    BC.States(i).Invariant.D = D;
    BC.States(i).Invariant.k = k;
    
    %iterate over transitions starting in State
    if isfield(State,'Trans')
        numTrans = length(State.Trans);
    else
        numTrans = 0;
    end
    for j = 1:numTrans
        Tran = State.Trans(j);
        
        %derive polytope for guard
        if convtype
            set = eq2set(Tran.guard.inequalities,Tran.guard.equalities, ...
                         BC.states,BC.outputsGlobal);
        else
            set = eq2set(Tran.guard.inequalities,Tran.guard.equalities,...
                         BC.states,[BC.outputsLocal BC.outputsGlobal]);
        end
        BC.States(i).Trans(j).guard.set = set;
        
        % derive linear representation of assignment
        [isLin,A,B,c,eqs,containsInput] = eq2linSys(Tran.reset.varNames,Tran.reset.expressions,...
                                                    BC.states,BC.inputs,'assignment');
        if(isLin)
            BC.States(i).Trans(j).reset.A = A;
            BC.States(i).Trans(j).reset.c = c;
            if containsInput
                % reset is input dependant
                BC.States(i).Trans(j).reset.B = B;
                BC.States(i).Trans(j).reset.hasInput = 1;
                BC.States(i).Trans(j).reset.inputDim = size(B,2);
            end
        else
            if strlength(eqs) ~= 0
                BC.States(i).Trans(j).reset.FormalEqs = eqs;
                % check if reset is input dependant
                if containsInput
                    BC.States(i).Trans(j).reset.hasInput = 1;
                    BC.States(i).Trans(j).reset.inputDim = length(BC.inputs);
                end
            else
                warning("Transition from ""%s"" to ""%s"" has a no reset function. Using identity...",...
                    State.name, BC.States(Tran.destination).name);
                BC.States(i).Trans(j).reset.A = diag(ones(size(BC.states)));
                BC.States(i).Trans(j).reset.c = zeros(size(BC.states));
            end
        end
    end
end

end