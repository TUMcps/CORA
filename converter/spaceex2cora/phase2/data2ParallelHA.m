function [functionName,HA] = data2ParallelHA(Data,path)
% data2ParallelHA - write text for hybrid automaton file
%
% Syntax:  
%    [functionName,HA] = data2ParallelHA(Data,path)
%
% Inputs:
%    Data - Automaton in structHA format
%    path - folder from which to generate auxillary files
%
% Outputs:
%    functionName - name of function
%    HA - text for hybridAutomaton instantiation
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

% Get Meta Information of the given Automaton
Components = Data.Components;
functionName = Data.name;
automaton_id = Data.componentID;

% Create main comments
functionStr = "function HA = " + functionName + "(~)";
dateComment = "%% Generated on " + datestr(date);
aCommentStr = padComment("Automaton created from Component '" + automaton_id + "'");
automatonStr = functionStr + newlines(3) + dateComment + newlines(2) + aCommentStr + newlines(2);

% Create Interface information comment
infoCommentStr = "%% Interface Specification:" + newline +...
    "%   This section clarifies the meaning of state, input & output dimensions" + newline +...
    "%   by showing their mapping to SpaceEx variable names. " + newlines(2);

% For each component list variable names of state & input space
numberOfComp = length(Components);
infoStr = "";
for comp = 1:numberOfComp
    Comp = Components{comp};
    compInfoStr = "% Component " + int2str(comp) + " (" + Comp.name + "):" + newline;
    
    % gather state names in string array
    stateNames = [Comp.states.name];
    % print to comment
    stateStr = "%  state x := [" + stateNames(1);
    if length(stateNames) > 1
        stateStr = stateStr + sprintf("; %s",stateNames(2:end));
    end
    stateStr = stateStr + "]" + newline;
    
    % gather input names in string array
    inputNames = [Comp.inputs.name];
    % print to comment
    inputStr = "%  input u := [" + inputNames(1);
    if length(inputNames) > 1
        inputStr = inputStr + sprintf("; %s",inputNames(2:end));
    end
    inputStr = inputStr + "]" + newline;

    % gather output names in string array
    if numberOfComp > 1
        % parallel
        outputStr = "%  output y := [";
        if ~isempty(Comp.outputsLocal)
            outputNames = [Comp.outputsLocal.name];
            % print to comment
            outputStr = outputStr + outputNames(1);
            if length(outputNames) > 1
                outputStr = outputStr + sprintf("; %s",outputNames(2:end));
            end
        end
        outputStr = outputStr + "]" + newlines(2);
    else
        % flat
        if ~isempty(Comp.outputsGlobal)
            % gather output names in string array
            outputNames = [Comp.outputsGlobal.name];
            % print to comment
            outputStr = "%  output y := [" + outputNames(1);
            if(length(outputNames) > 1)
                outputStr = outputStr + sprintf("; %s",outputNames(2:end));
            end
            outputStr = outputStr + "]" + newlines(2);
        else
            outputStr = newline;
        end
    end
    
    infoStr = infoStr + compInfoStr + stateStr + inputStr + outputStr;
end

automatonStr = automatonStr + infoCommentStr + infoStr;

componentStr = "";
% for stateBinds of pHA
totalStates = 0;

% For each component in automaton
for comp = 1:numberOfComp
    Comp = Components{comp};
    
    %Get Meta Information for the Component "comp"
    component_id = Comp.name;
    States = Comp.States;
    
    % Write Comment for Component "comp"
    cCommentStr = padComment("Component " + component_id) + newlines(2);
    % Append it to component String
    componentStr = componentStr + cCommentStr;

    % re-initialize 'loc' variable, append to component string
    if comp > 1
        componentStr = componentStr + "loc = {};" + newlines(2);
    end
    
    
    % For each state in component
    numberOfStates = length(States);
    for state = 1:numberOfStates
        State = States(state);
        
        % Write Comment for State "state"
        sCommentStr = padComment("State " + State.name) + newlines(2);
        stateStr = sCommentStr;
        
        % Give original equation as comment
        dynamicsC = text2comment("flow equation:" + newline + State.Flow.Text) + newline;
        if isfield(State.Flow,'A')
            % Get information for linear system
            linSysA = printMatrixConverter(State.Flow.A);
            linSysAStr = "dynA = ..." + newline + linSysA + ";" + newline;
            linSysB = printMatrixConverter(State.Flow.B);
            linSysBStr = "dynB = ..." + newline + linSysB + ";" + newline;
            linSysc = printMatrixConverter(State.Flow.c);
            linSyscStr = "dync = ..." + newline + linSysc + ";" + newline;
            % Get information about outputs from Invariant of linear system
            if (numberOfComp == 1 && isempty(Comp.outputsGlobal)) || ...
                (numberOfComp > 1 && isempty(Comp.outputsLocal))
                % no outputs, no equation of the form y = Cx + Du + k
                dynamicsStr = dynamicsC + linSysAStr + linSysBStr + linSyscStr + ...
                    "dynamics = linearSys(dynA, dynB, dync);" + newlines(2);
            else
                % include output equation y = Cx + Du + k
                linSysC = printMatrixConverter(State.Flow.C);
                linSysCStr = "dynC = ..." + newline + linSysC + ";" + newline;
                linSysD = printMatrixConverter(State.Flow.D);
                linSysDStr = "dynD = ..." + newline + linSysD + ";" + newline;
                linSysk = printMatrixConverter(State.Flow.k);
                linSyskStr = "dynk = ..." + newline + linSysk + ";" + newline;
                dynamicsStr = dynamicsC + linSysAStr + linSysBStr + linSyscStr + ...
                    linSysCStr + linSysDStr + linSyskStr + ...
                    "dynamics = linearSys(dynA, dynB, dync, dynC, dynD, dynk);" + newlines(2);
            end
        else
            % choose name for dynamics function
            if numberOfComp==1
                % simplify names for monolithic automata
                nonlinName = sprintf("%s_Loc%d_FlowEq",functionName,state);
            else
                nonlinName = sprintf("%s_Comp%d_Loc%d_FlowEq",functionName,comp,state);
            end
            
            % find dynamics of system
            statedims = num2str(length(Comp.states));
            inputdims = num2str(length(Comp.inputs));
            
            printDynamicsFile(path,nonlinName,State.Flow.FormalEqs,"flow");
            
            dynamicsStr = dynamicsC + "dynamics = nonlinearSys(@" + ...
                          nonlinName + "," + statedims + "," + ...
                          inputdims + "); " + newlines(2);
        end
        
        % Get information for Invariant
        InvText = State.Invariant.Text;
        if isnumeric(State.Invariant.set) && isempty(State.Invariant.set)
            str1 = "";
            str2 = "[];";
            InvText = "no invariant given";
        elseif isa(State.Invariant.set,'mptPolytope')
            [str1,str2] = mptPolytopeString(State.Invariant.set);
        elseif isa(State.Invariant.set,'levelSet')
            [str1,str2] = levelSetString(State.Invariant.set);
        else
            throw(CORAerror('CORA:converterIssue',...
                'Invariant has to be either empty, an mptPolytope, or a levelSet'));
        end

        % Write String for Invariant
        invariantC = text2comment("invariant equation:" + newline + InvText) + newline;
        invariantStr = invariantC + str1 + "inv = " + str2 + newlines(2);
        
        transitionStr = "trans = {};" + newline;
        % For each Transition
        Trans = State.Trans;
        numberOfTrans = length(Trans);
        for trans = 1:numberOfTrans
            Tran = Trans(trans);
            
            % Get information for destination for Transition "trans"
            transDestination = num2str(Tran.destination);
            
            % Get Information for Reset for Transition "trans"
            if isfield(Tran.reset,'A')
                % linear reset
                resetA = printMatrixConverter(Tran.reset.A);
                resetAStr = "resetA = ..." + newline + resetA + ";" + newline;
                resetc = printMatrixConverter(Tran.reset.c);
                resetcStr = "resetc = ..." + newline + resetc + ";" + newline;
                if isfield(Tran.reset,'B')
                    resetB = printMatrixConverter(Tran.reset.B);
                    resetBStr = "resetB = ..." + newline + resetB + ";" + newline;
                end
                
                % Write Reset String
                tranResetText = Tran.reset.Text;
                if tranResetText == ""
                    tranResetText = "no reset equation given";
                end
                resetComment = text2comment("reset equation:" + newline + tranResetText) + newline;
                if isfield(Tran.reset,'B')
                    resetStr = resetComment + resetAStr + resetBStr + resetcStr + ...
                        "reset = struct('A', resetA, 'B', resetB, 'c', resetc);" + newlines(2);
                else
                    resetStr = resetComment + resetAStr + resetcStr + ...
                        "reset = struct('A', resetA, 'c', resetc);" + newlines(2);
                end
            else
                % nonlinear reset
                               
                % choose name of function
                if numberOfComp == 1
                    resetFuncName = sprintf("%s_Loc%d_Trans%d_ResetEq%d",...
                        functionName,state,trans);
                else
                    resetFuncName = sprintf("%s_Comp%d_Loc%d_Trans%d_ResetEq%d",...
                        functionName,comp,state,trans);
                end

                % generate function file for nonlinear reset
                printDynamicsFile(path,resetFuncName,Tran.reset.FormalEqs,"reset");
                tranResetText = Tran.reset.Text;
                if tranResetText == ""
                    tranResetText = "no reset equation given";
                end
                
                % Write Reset String
                resetComment = text2comment("reset equation:" + ...
                    newline + tranResetText) + newline;
                resetStr = resetComment + "reset = struct('f', @" + ...
                        resetFuncName + ");" + newlines(2);
            end
            
            % Get Information for Guards for Transition "trans"
            tranGuardText = Tran.guard.Text;
            if tranGuardText == ""
                tranGuardText = "no guard set given";
            end
            if isnumeric(Tran.guard.set) && isempty(Tran.guard.set)
                % empty guard set -> immediate transition
                str1 = "";
                str2 = "[];";

            elseif isa(Tran.guard.set,'mptPolytope')
                % intersect guard with invariant
                try
                    G = State.Invariant.set & Tran.guard.set;                    
                catch
                    G = Tran.guard.set;
                end

                % check if guard can be represented as hyperplane
                res = false;
                if isa(G,'mptPolytope')
                    G = removeRedundancies(G,'all');
                    [res,ch] = isConHyperplane(G);
                end
                
                if isa(G,'levelSet')
                    [str1,str2] = levelSetString(G);
                elseif isa(G,'mptPolytope')
                    if res
                        [str1,str2] = conHyperplaneString(ch);
                    else
                        [str1,str2] = mptPolytopeString(G);
                    end
                else
                    throw(CORAerror('CORA:converterIssue',...
                        'Unexpected set representation in conversion of guard set.'));
                end
            elseif isa(Tran.guard.set,'levelSet')
                % intersect guard with invariant
                try
                    G = State.Invariant.set & Tran.guard.set;                    
                catch
                    G = Tran.guard.set;
                end
                [str1,str2] = levelSetString(G);
            else
                throw(CORAerror('CORA:converterIssue',...
                    ['Guard set has to be empty, a conHyperplane, '...
                    'an mptPolytope or a levelSet.'])); 
            end
            
            % Write Guard String
            guardC = text2comment("guard equation:" + newline + tranGuardText) + newline;
            guardStr = guardC + str1 + "guard = " + str2 + newlines(2);
            
            
            % Write Transition String, include only non-empty labels
            if strlength(Tran.label) > 0
                transStr = "trans{" + num2str(trans) + ...
                    "} = transition(guard, reset, " + ...
                    transDestination + ", '" + Tran.label + "');" + newlines(2);
            else
                transStr = "trans{" + num2str(trans) + ...
                    "} = transition(guard, reset, " + ...
                    transDestination + ");" + newlines(2);
            end
            % Append Transition string
            transitionStr = transitionStr + resetStr + guardStr + transStr;
            
        end
        
        % Write State String
        if States(state).name == ""
            locName = "S" + state;
        else
            locName = States(state).name;
        end
        locStr = "loc{" + num2str(state) + "} = location('" + locName + ...
            "', inv, trans, dynamics);" + newlines(4);
        % Append State String
        stateStr = stateStr + dynamicsStr + invariantStr + transitionStr + locStr;
        % Append State String to Component String
        componentStr = componentStr + stateStr;
        
    end
    
    if numberOfComp > 1
        % only parallel hybrid automata here...
        parallelCompStr = text2comment("composition: hybrid automaton and input binds") + ...
            newline + "comp{" + num2str(comp) + "} = hybridAutomaton(loc);" + newlines(2);
        
        % declare inputBinds
        % default string for comp without inputs
        inputBindStr = "iBinds{" + comp + "} = [];" + newlines(2);
        % syntax of inputBinds: mx2 array, where
        %   m - number of input arguments to current component
        % with (m,1) being either 0 (global) or the component of origin
        % and (m,2) being the number of the global input (if (m,1) = 0)
        %           or the number of output of the component of origin

        if strcmp(Comp.inputs(1).name,'uDummy')
            % only dummy input generated by conversion... use first global
            % input as input bind (no effect as all input matrices are 0)
            inputBindStr = "% only dummy input" + newline + ...
                "iBinds{" + comp + "} = [0 1];" + newlines(2);

        else
            % actual (meaningful) inputs given

            inputBindStr = "% input names: " + strjoin([Comp.inputs.name],", ") + newline + ...
                "iBinds{" + comp + "} = [[";

            % loop over all inputs
            for inp=1:length(Comp.inputs)
                if inp > 1
                    inputBindStr = inputBindStr + ";[";
                end
                inputName = Comp.inputs(inp).name;
                % default: input comes from composed system ('global')
                origin = 0;
                outputOfOrigin = 1;
                % search for inputName in all other components
                for c=1:length(Components)
                    idx = strcmp([Components{c}.outputsLocal.name],inputName);
                    if any(idx)
                        origin = c;
                        % name should only occur only once
                        if nnz(idx) > 1
                            throw(CORAerror('CORA:converterIssue',...
                                ['The input name of an input bind could not be resolved '...
                                'because the name occurs twice in another component.']));
                        end
                        % index in list of outputs of other component
                        outputOfOrigin = find(idx,1,'first');
                        % continue with next row in input binds
                        break
                    end
                end
                % write bind
                inputBindStr = inputBindStr + num2str(origin) + "," +...
                    num2str(outputOfOrigin) + "]";
            end
            inputBindStr = inputBindStr + "];" + newlines(2);
        end
        
        % put all together
        componentStr = componentStr + parallelCompStr + inputBindStr;
    end
    
end


if numberOfComp == 1
    % If the number of Components is 1, we have a flat automaton   
    aStr = "HA = hybridAutomaton(loc);" + newlines(2);
else
    % If the number of Components is > 1, we have a parallel hybrid automaton
    aStr = "HA = parallelHybridAutomaton(comp,iBinds);" + newlines(2);
    
    if all(strcmp(cellfun(@(cpnt) cpnt.inputs.name, Components), 'uDummy'))
    % if no inputs in system, CORA needs an input
    % generate a global input on the first component
    % this will have no effect as the input matrices should all be 0
        componentStr = componentStr + ...
            "% no inputs given, global input for CORA, no effect as B = 0" + ...
            newline + "iBinds{1} = [0,1];" + newlines(2);
    end
end

%optionStr = padComment("Options");

HA = automatonStr + componentStr + aStr + newline + "end";

end


% Auxiliary Functions -----------------------------------------------------

function str = padComment(comment,maxLineLength)
%pads comment left & right with dashes to desired length and prefixes "%"

    if nargin < 2
        maxLineLength = 75;
    end
    
    lenComment = strlength(comment);
    lenLeft = floor((maxLineLength - lenComment)/2) - 1;
    lenRight = maxLineLength - lenLeft - lenComment;
    
    str = "%" + repmat('-',1,lenLeft-1) + comment + repmat('-',1,lenRight);

end

function str = newlines(lines)
% fast way to write newline() + newline() + ...
    str = string(repmat(newline(),1,lines));
end

function str = text2comment(text)
% transform possibly multi-line text to comment
% format in:
%   "line1
%    line2
%    line3"
% format out:
%   "%% line1
%    %   line2
%    %   line3"
    str = "%% " + strrep(text,newline,newline + "%   ");
end

function [str1,str2] = mptPolytopeString(set)
% generates the string that constructs the mptPolytope

    A = printMatrixConverter(set.P.A);
    AStr = "A = ..." + newline + A + ";" + newline;
    b = printMatrixConverter(set.P.b);
    bStr = "b = ..." + newline + b + ";" + newline;
    OptStr = "polyOpt = struct('A', A, 'b', b";
    if ~isempty(set.P.Ae)
        % if invariant is a polyhedron, add additional parameters
        Ae = printMatrixConverter(set.P.Ae);
        AeStr = "Ae = ..." + newline + Ae + ";" + newline;
        be = printMatrixConverter(set.P.be);
        beStr = "be = ..." + newline + be + ";" + newline;
        OptStr = OptStr + ",'Ae', Ae, 'be', be";
    else
        AeStr = "";
        beStr = "";
    end
    str1 = AStr + bStr + AeStr + beStr + ...
                OptStr + ");" + newline;

    str2 = 'mptPolytope(polyOpt);';
end

function [str1,str2] = levelSetString(set)
% generates the string that constructs the levelSet
    
    % generate string for constructing the variable vector
    varStr = "vars = sym('x',[" + num2str(length(set.vars)) + ",1]);" + newline;
    
    varStr = varStr + "syms ";
    for i = 1:length(set.vars)
        if i ~= length(set.vars)
            varStr = varStr + "x" + num2str(i) + " ";
        else
            varStr = varStr + "x" + num2str(i) + ";";
        end
    end
    
    % generate string for constructing the symbolic equations
    x = sym('x',[length(set.vars),1]);
    eq = set.funHan(x);
    
    if length(eq) == 1
        eqStr = "eq = " + string(eq) + ";";
    else
        eqStr = "eq = [";
        for i = 1:length(eq)
            if i ~= length(eq)
                eqStr = eqStr + string(eq(i)) + "; ..." + newline;
            else
                eqStr = eqStr + string(eq(i));
            end
        end
        eqStr = eqStr + "];";
    end
    
    % generate string for the comparison operators
    if ~iscell(set.compOp)
        compStr = "compOp = '" + set.compOp + "';"; 
    else
        compStr = "compOp = {";
        for i = 1:length(set.compOp)
           if i ~= length(set.compOp)
                compStr = compStr + "'" + set.compOp{i} + "',";
           else
                compStr = compStr + "'" + set.compOp{i} + "'};";
           end
        end
    end
    
    % generate overall string
    str1 = varStr + newline + eqStr + newline + compStr + newlines(2);
    str2 = "levelSet(eq,vars,compOp);";

end

function [str1,str2] = conHyperplaneString(set)
% generates the string that constructs the conHyperplane object

    % hyperplane equation c*x = d
    c = printMatrixConverter(set.h.c);
    cStr = "c = " + c + ";" + newline;
    dStr = "d = " + num2str(set.h.d) + ";";
    
    % inequality constraints C*x <= D
    if ~isempty(set.C)
        C = printMatrixConverter(set.C);
        Cstr = "C = ..." + newline + C + ";" + newline;
        D = printMatrixConverter(set.d);
        Dstr = "D = " + D + ";";
        
        str1 = cStr + dStr + Cstr + Dstr + newlines(2);
        str2 = "conHyperplane(c,d,C,D);";
    else
        str1 = cStr + dStr + newlines(2);
        str2 = "conHyperplane(c,d);";
    end
end

%------------- END OF CODE --------------
