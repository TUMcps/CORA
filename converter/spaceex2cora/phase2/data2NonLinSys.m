function [functionName,sys] = data2NonLinSys(Data, path)
% data2NonLinSys - converts Data into a linearSys or nonlinearSys object
%
% Syntax:  
%    [functionName,sys] = data2NonLinSys(Data, path)
%
% Input:
%    Data: Automaton in structHA format
%    path: folder from which to generate auxiliary files
% Output:
%     depending on system dynamics:
%       linearSys object or nonlinearSys object
%
% Other m-files required: none
% Subfunctions: ADD IF NECESSARY
% MAT-files required: none
%
% See also: none

% Author: Mark Wetzlinger
% Written: 11-January-2019
% Last update: --- 
% Last revision: ---

%------------- BEGIN CODE --------------


% Get Meta Information of the given Automaton
Comp = Data.Components{1,1};
functionName = Data.name;
automaton_id = Data.componentID;

% Create main comments
functionStr = "function sys = " + functionName + "(~)";
dateComment = "%% Generated on " + datestr(date);
aCommentStr = padComment("Automaton created from Component '" + automaton_id + "'");
automatonStr = functionStr + newlines(3) + dateComment + newlines(2) + aCommentStr + newlines(2);

% Create Interface information comment
infoCommentStr = "%% Interface Specification:" + newline +...
    "%   This section clarifies the meaning of state, input & output dimensions" + newline +...
    "%   by showing their mapping to SpaceEx variable names. " + newlines(2);

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
    outputStr = outputStr + "]" + newlines(2);
else
    outputStr = newline;
end

infoStr = infoStr + compInfoStr + stateStr + inputStr + outputStr;

automatonStr = automatonStr + infoCommentStr + infoStr;

%Get Meta Information for the Component "comp"
component_id = Comp.name;
State = Comp.States(1,1);
    
% Write Comment for Component "comp"
cCommentStr = padComment("Component " + component_id) + newlines(2);
% Append it to component String
componentStr = cCommentStr;

% Write Comment for State "state"
sCommentStr = padComment("State " + State.name) + newlines(2);
stateStr = sCommentStr;

% Give original equation as comment
dynamicsC = text2comment("equation:" + newline + State.Flow.Text) + newline;
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
            "sys = linearSys(dynA, dynB, dync);" + newlines(2);
    else
        % include output equation y = Cx + Du + k
        linSysC = printMatrixConverter(State.Invariant.C);
        linSysCStr = "dynC = ..." + newline + linSysC + ";" + newline;
        linSysD = printMatrixConverter(State.Invariant.D);
        linSysDStr = "dynD = ..." + newline + linSysD + ";" + newline;
        linSysk = printMatrixConverter(State.Invariant.k);
        linSyskStr = "dynk = ..." + newline + linSysk + ";" + newline;
        dynamicsStr = dynamicsC + linSysAStr + linSysBStr + linSyscStr + ...
            linSysCStr + linSysDStr + linSyskStr + newlines(2) + ...
            "sys = linearSys(dynA, dynB, dync, dynC, dynD, dynk);" + newlines(2);
    end
else
    % choose name for dynamics function
    % simplify names for monolithic automata
    nonlinName = sprintf("%s_St%d_FlowEq",functionName,1);

    % find dynamics of system
    statedims = num2str(length(Comp.states));
    inputdims = num2str(length(Comp.inputs));

    printDynamicsFile(path,nonlinName,State.Flow.FormalEqs);

    dynamicsStr = dynamicsC + "sys = nonlinearSys(@" + nonlinName + ...
                  "," + statedims + "," + inputdims +"); " + newlines(2);
end

% Append State String to Component String
componentStr = componentStr + stateStr + dynamicsStr;

%optionStr = padComment("Options");


sys = automatonStr + componentStr + newline + "end";

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

%------------- END OF CODE -------------