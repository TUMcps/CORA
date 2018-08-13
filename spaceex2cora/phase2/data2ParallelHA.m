function [functionName,HA] = data2ParallelHA( Data, path)

%INPUT :
%    Data : Automaton in structHA formate
%    path: folder from which to generate auxillary files
%OUTPUT :
%     HA: Automaton in CORA format as a string
%


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
    "%   This section clarifies the meaning of state & input dimensions" + newline+...
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
    if(length(stateNames) > 1)
        stateStr = stateStr + sprintf("; %s",stateNames(2:end));
    end
    stateStr = stateStr + "]" + newline;
    
    % gather input names in string array
    inputNames = [Comp.inputs.name];
    % print to comment
    inputStr = "%  input u := [" + inputNames(1);
    if(length(inputNames) > 1)
        inputStr = inputStr + sprintf("; %s",inputNames(2:end));
    end
    inputStr = inputStr + "]" + newlines(2);
    
    infoStr = infoStr + compInfoStr + stateStr + inputStr;
end

automatonStr = automatonStr + infoCommentStr + infoStr;

% For each component in automaton
for comp = 1:numberOfComp
    
    Comp = Components{comp};
    
    %Get Meta Information for the Component "comp"
    component_id = Comp.name;
    States = Comp.States;
    
    % Write Comment for Component "comp"
    cCommentStr = padComment("Component " + component_id) + newlines(2);
    % Append it to component String
    componentStr = cCommentStr;
    
    
    % For each state in component
    numberOfStates = length(States);
    for state = 1:numberOfStates
        State = States(state);
        
        % Write Comment for State "state"
        sCommentStr = padComment("State " + State.name) + newlines(2);
        stateStr = sCommentStr;
        
        % Give original equation as comment
        dynamicsC = text2comment("equation:" + newline + State.Flow.Text) + newline;
        if isfield(State.Flow,'A')
            % Get information for linear System
            linSysA = printMatrixConverter(State.Flow.A);
            linSysAStr = "dynA = ..." + newline + linSysA + ";" + newline;
            linSysB = printMatrixConverter(State.Flow.B);
            linSysBStr = "dynB = ..." + newline + linSysB + ";" + newline;
            linSysc = printMatrixConverter(State.Flow.c);
            linSyscStr = "dync = ..." + newline + linSysc + ";" + newline;
            dynamicsStr = dynamicsC + linSysAStr + linSysBStr + linSyscStr + ...
                    "dynamics = linearSys('linearSys', dynA, dynB, dync);" + newlines(2);
        else
            % choose name for dynamics function
            if numberOfComp==1
                %siplify names for monolithic automata
                nonlinName = sprintf("%s_St%d_FlowEq",functionName,state);
            else
                nonlinName = sprintf("%s_C%d_St%d_FlowEq",functionName,comp,state);
            end
            
            dynOptStr = "dynOpt = struct('tensorOrder',1);" + newline;
            
            %find dynamics of system
            statedims = num2str(length(Comp.states));
            inputdims = num2str(length(Comp.inputs));
            
            printDynamicsFile(path,nonlinName,State.Flow.FormalEqs);
            
            dynamicsStr = dynamicsC + dynOptStr + "dynamics = nonlinearSys(" +...
                statedims + "," + inputdims + ",@" + nonlinName + ",dynOpt); " + newlines(2);
        end
        
        % Get information for Invariant
        invariantA = printMatrixConverter(State.Invariant.A);
        invariantAStr = "invA = ..." + newline + invariantA + ";" + newline;
        invariantb = printMatrixConverter(State.Invariant.b);
        invariantbStr = "invb = ..." + newline + invariantb + ";" + newline;
        invOptStr = "invOpt = struct('A', invA, 'b', invb";
        if ~isempty(State.Invariant.Ae)
            % if invariant is a polyhedron, add additional parameters
            invariantAe = printMatrixConverter(State.Invariant.Ae);
            invariantAeStr = "invAe = ..." + newline + invariantAe + ";" + newline;
            invariantbe = printMatrixConverter(State.Invariant.be);
            invariantbeStr = "invbe = ..." + newline + invariantbe + ";" + newline;
            invOptStr = invOptStr + ",'Ae', invAe, 'be', invbe";
        else
            invariantAeStr = "";
            invariantbeStr = "";
        end
        invOptStr = invariantAStr + invariantbStr + invariantAeStr + invariantbeStr + ...
                    invOptStr + ");" + newline;
        
        % Write String for Invariant
        invariantC = text2comment("equation:" + newline + State.Invariant.Text) + newline;
        invariantStr = invariantC + invOptStr + "inv = mptPolytope(invOpt);" + newlines(2);
        
        transitionStr = "trans = {};" + newline;
        % For each Transition
        Trans = State.Trans;
        numberOfTrans = length(Trans);
        for trans = 1:numberOfTrans
            Tran = Trans(trans);
            
            % Get information for destination for Transition "trans"
            transDestination = num2str(Tran.destination);
            
            % Get Information for Reset for Transition "trans"
            resetA = printMatrixConverter(Tran.reset.A);
            resetAStr = "resetA = ..." + newline + resetA + ";" + newline;
            resetb = printMatrixConverter(Tran.reset.b);
            resetbStr = "resetb = ..." + newline + resetb + ";" + newline;
            
            % Write Reset String
            resetC = text2comment("equation:" + newline + Tran.reset.Text) + newline;
            resetStr = resetC + resetAStr + resetbStr + ...
                    "reset = struct('A', resetA, 'b', resetb);" + newlines(2);
            
            % Get Information for Guards for Transition "trans"
            guardA = printMatrixConverter(Tran.guard.A);
            guardAStr = "guardA = ..." + newline + guardA + ";" + newline;
            guardb = printMatrixConverter(Tran.guard.b);
            guardbStr = "guardb = ..." + newline + guardb + ";" + newline;
            guardOptStr = "guardOpt = struct('A', guardA, 'b', guardb";
            if ~isempty(Tran.guard.Ae)
                % if guard is a polyhedron, add additional parameters
                guardAe = printMatrixConverter(Tran.guard.Ae);
                guardAeStr = "guardAe = ..." + newline + guardAe + ";" + newline;
                guardbe = printMatrixConverter(Tran.guard.be);
                guardbeStr = "guardbe = ..." + newline + guardbe + ";" + newline;
                guardOptStr = guardOptStr + ", 'Ae', guardAe, 'be', guardbe";
            else
                guardAeStr = "";
                guardbeStr = "";
            end
            guardOptStr = guardAStr + guardbStr + guardAeStr + guardbeStr + ...
                        guardOptStr + ");" + newline;
            
            % Write Guard String
            guardC = text2comment("equation:" + newline + Tran.guard.Text) + newline;
            guardStr = guardC + guardOptStr + "guard = mptPolytope(guardOpt);" + newlines(2);
            
            % Write Transition String
            transStr = "trans{" + num2str(trans) + "} = transition(guard, reset, " +...
                        transDestination + ", 'dummy', 'names');" + newlines(2);
            % Append Transition string
            transitionStr = transitionStr + resetStr + guardStr + transStr;
            
        end
        
        % Write State String
        locStr = "loc{" + num2str(state) + "} = location('S" + num2str(state) +...
                    "'," + num2str(state) + ", inv, trans, dynamics);" + newlines(4);
        % Append State String
        stateStr = stateStr + dynamicsStr + invariantStr + transitionStr + locStr;
        % Append State String to Component String
        componentStr = componentStr + stateStr;
        
    end
    
end
if numberOfComp >1
% If the number of Components is > 1, we have a parallel hybrid Automaton
    aStr = "HA = parallelHybridAutomaton(comp);" + newlines(2); 
else
% If the number of Components is 1, we have a monolithic Automaton   
    aStr = "HA = hybridAutomaton(loc);" + newlines(2);
end

%optionStr = padComment("Options");


HA = automatonStr + componentStr + aStr + newline + "end";

end

%-----------STRING HELPER FUNCTIONS-------------

function str = padComment(comment,maxLineLength)
%pads comment left & right with dashes to desired length and prefixes "%"

if(nargin<2)
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

% transform possibly multi-line text to comment
function str = text2comment(text)
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
