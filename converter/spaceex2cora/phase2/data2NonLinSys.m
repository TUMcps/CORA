function sys = data2NonLinSys(data,functionName,path)
% data2NonLinSys - converts Data into a linearSys or nonlinearSys object
%
% Syntax:
%    sys = data2NonLinSys(data,functionName,path)
%
% Inputs:
%    data - linear/nonlinear system in structHA format
%    functionName - name of CORA model file
%    path - folder from which to generate auxiliary files
%
% Outputs:
%    sys - linearSys object or nonlinearSys object (depending on dynamics)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       11-January-2019
% Last update:   --- 
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% Get Meta Information of the given Automaton
Comp = data.components{1,1};
automaton_id = data.componentID;

% Create main comments
functionStr = "function sys = " + functionName + "(~)";
dateComment = "%% Generated on " + datestr(date);
aCommentStr = aux_padComment("Automaton created from Component '" + automaton_id + "'");
automatonStr = functionStr + aux_newlines(3) + dateComment + aux_newlines(2) + aCommentStr + aux_newlines(2);

% Create Interface information comment
infoCommentStr = "%% Interface Specification:" + newline +...
    "%   This section clarifies the meaning of state, input & output dimensions" + newline +...
    "%   by showing their mapping to SpaceEx variable names. " + aux_newlines(2);

% For each component list variable names of state & input space
infoStr = "";
compInfoStr = "% Component (" + Comp.name + "):" + newline;
    
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
inputStr = inputStr + "]" + newline;
    
if ~isempty(Comp.outputsGlobal)
    % gather output names in string array
    outputNames = [Comp.outputsGlobal.name];
    % print to comment
    outputStr = "%  output y := [" + outputNames(1);
    if(length(outputNames) > 1)
        outputStr = outputStr + sprintf("; %s",outputNames(2:end));
    end
    outputStr = outputStr + "]" + aux_newlines(2);
else
    outputStr = newline;
end

infoStr = infoStr + compInfoStr + stateStr + inputStr + outputStr;

automatonStr = automatonStr + infoCommentStr + infoStr;

%Get Meta Information for the Component "comp"
component_id = Comp.name;
State = Comp.States(1,1);
    
% Write Comment for Component "comp"
cCommentStr = aux_padComment("Component " + component_id) + aux_newlines(2);
% Append it to component String
componentStr = cCommentStr;

% Write Comment for State "state"
sCommentStr = aux_padComment("State " + State.name) + aux_newlines(2);
stateStr = sCommentStr;

% Give original equation as comment
dynamicsC = aux_text2comment("equation:" + newline + State.Flow.Text) + newline;
if isfield(State.Flow,'A')
    % Get information for linear system
    linSysA = printMatrixConverter(State.Flow.A);
    linSysAStr = "dynA = ..." + newline + linSysA + ";" + newline;
    linSysB = printMatrixConverter(State.Flow.B);
    linSysBStr = "dynB = ..." + newline + linSysB + ";" + newline;
    linSysc = printMatrixConverter(State.Flow.c);
    linSyscStr = "dync = ..." + newline + linSysc + ";" + newline;
    % Get information about outputs from Invariant of linear system
    if isempty(Comp.outputsGlobal)
        % no outputs, no equation of the form y = Cx + Du + k
        dynamicsStr = dynamicsC + linSysAStr + linSysBStr + linSyscStr + ...
            "sys = linearSys(dynA, dynB, dync);" + aux_newlines(2);
    else
        % include output equation y = Cx + Du + k
        linSysC = printMatrixConverter(State.Flow.C);
        linSysCStr = "dynC = ..." + newline + linSysC + ";" + newline;
        linSysD = printMatrixConverter(State.Flow.D);
        linSysDStr = "dynD = ..." + newline + linSysD + ";" + newline;
        linSysk = printMatrixConverter(State.Flow.k);
        linSyskStr = "dynk = ..." + newline + linSysk + ";" + newline;
        dynamicsStr = dynamicsC + linSysAStr + linSysBStr + linSyscStr + ...
            linSysCStr + linSysDStr + linSyskStr + aux_newlines(2) + ...
            "sys = linearSys(dynA, dynB, dync, dynC, dynD, dynk);" + aux_newlines(2);
    end
else
    % choose name for dynamics function
    % simplify names for monolithic automata
    nonlinName = sprintf("%s_St%d_FlowEq",functionName,1);

    % find dynamics of system
    statedims = num2str(length(Comp.states));
    inputdims = num2str(length(Comp.inputs));

    printDynamicsFile(path,nonlinName,State.Flow.FormalEqs,'flow');

    dynamicsStr = dynamicsC + "sys = nonlinearSys(@" + nonlinName + ...
                  "," + statedims + "," + inputdims +"); " + aux_newlines(2);
end

% Append State String to Component String
componentStr = componentStr + stateStr + dynamicsStr;

%optionStr = aux_padComment("Options");

sys = automatonStr + componentStr + newline + "end";

end


% Auxiliary functions -----------------------------------------------------

function str = aux_padComment(comment,maxLineLength)
% pads comment left & right with dashes to desired length and prefixes "%"
    if(nargin<2)
        maxLineLength = 75;
    end
    
    lenComment = strlength(comment);
    lenLeft = floor((maxLineLength - lenComment)/2) - 1;
    lenRight = maxLineLength - lenLeft - lenComment;
    
    str = "%" + repmat('-',1,lenLeft-1) + comment + repmat('-',1,lenRight);
end

function str = aux_newlines(lines)
    % fast way to write newline() + newline() + ...
    str = string(repmat(newline(),1,lines));
end

function str = aux_text2comment(text)
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

% ------------------------------ END OF CODE ------------------------------
