function HA = data2ParallelHA(data,functionName,resultpath)
% data2ParallelHA - write text for hybrid automaton file
%
% Syntax:
%    HA = data2ParallelHA(data,functionName,resultpath)
%
% Inputs:
%    data - automaton in structHA format
%    functionName - name of CORA model file
%    resultpath - folder from which to generate auxiliary files
%
% Outputs:
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

% Authors:       ???, Mark Wetzlinger
% Written:       ---
% Last update:   ---
% Last revision: 11-January-2023 (MW, restructure code)

% ------------------------------ BEGIN CODE -------------------------------

% get meta information of the given automaton
components = data.components;
automaton_id = data.componentID;

% create main comments
functionStr = "function HA = " + functionName + "(~)";
dateComment = "%% Generated on " + datestr(datetime);
createdStr = aux_padComment("Automaton created from Component '" + automaton_id + "'");
automatonStr = functionStr + aux_newlines(3) + dateComment + aux_newlines(2) + createdStr + aux_newlines(2);

% create interface information comment
infoCommentStr = "%% Interface Specification:" + newline +...
    "%   This section clarifies the meaning of state, input & output dimensions" + newline +...
    "%   by showing their mapping to SpaceEx variable names. " + aux_newlines(2);

% number of components
nrComp = length(components);
% resulting automaton: flat or parallel?
isFlatHA = nrComp == 1;


% init information string (names of states, inputs, outputs for every
% component of the resulting automaton)
infoStr = aux_infoStr(components);

% expand full string
automatonStr = automatonStr + infoCommentStr + infoStr;

% init component string
componentStr = "";

% loop over each component in automaton
for iComp = 1:nrComp
    % read out i-th component
    comp = components{iComp};
    
    % write comment for component "comp"
    cCommentStr = aux_padComment("Component " + comp.name) + aux_newlines(2);
    % Append it to component String
    componentStr = componentStr + cCommentStr;

    % initialize 'loc' variable, (technically, we only need this for
    % iComp > 1, but add it everywhere for consistency)
    componentStr = componentStr + "clear loc" + aux_newlines(2);
    
    % shorter access to locations of i-th component
    locs = comp.States;
    
    % loop over each location
    nrLocs = length(locs);
    for iLoc = 1:nrLocs
        % shorter access to i-th location
        loc = locs(iLoc);
        % loop over all outgoing transitions
        allTrans = loc.Trans;
        
        % general information about i-th location
        locStr = aux_padComment("Location " + loc.name) + aux_newlines(2);
        
        % describe flow dynamics
        dynamicsStr = aux_dynamicsStr(loc,iLoc,comp,iComp,isFlatHA,resultpath,functionName);
        
        % describe invariant
        invariantStr = aux_invariantStr(loc,allTrans,isFlatHA);
        
        % init transition string
        transitionStr = "trans = transition();" + newline;

        for iTrans = 1:length(allTrans)
            % shorter access to i-th transition
            trans = allTrans(iTrans);
            
            % read out target location
            targetStr = num2str(trans.destination);

            % describe reset equation
            resetStr = aux_resetStr(trans,iComp,iLoc,iTrans,isFlatHA,resultpath,functionName);
            
            % describe guard set
            guardStr = aux_guardStr(trans,loc);

            % concatenate transition string incl. synchronization labels
            % (only for conversion to parallel hybrid automaton)
            if ~isFlatHA && strlength(trans.label) > 0
                transStr = "trans(" + num2str(iTrans) + ...
                    ") = transition(guard, reset, " + ...
                    targetStr + ", '" + trans.label + "');" + aux_newlines(2);
            else
                transStr = "trans(" + num2str(iTrans) + ...
                    ") = transition(guard, reset, " + ...
                    targetStr + ");" + aux_newlines(2);
            end
            % full transition string
            transitionStr = transitionStr + resetStr + guardStr + transStr;
            
        end
        
        % name of location
        if locs(iLoc).name == ""
            locName = "S" + iLoc;
        else
            locName = locs(iLoc).name;
        end
        locNameStr = "loc(" + num2str(iLoc) + ") = location('" + locName + ...
            "', inv, trans, dynamics);" + aux_newlines(4);

        % append location/dynamics/invariant/transition/location name
        locStr = locStr + dynamicsStr + invariantStr + transitionStr + locNameStr;
        % append to component string
        componentStr = componentStr + locStr;
        
    end
    
    % append input binds string (only for parallel hybrid automata)
    if ~isFlatHA
        inputBindsStr = aux_inputBindsStr(comp,iComp,components);
        componentStr = componentStr + inputBindsStr;
    end
    
end

% string to instantiate flat/parallel automaton
if isFlatHA
    % instantiate hybridAutomaton object
    objectStr = "HA = hybridAutomaton(loc);";
else
    % instantiate parallelHybridAutomaton object
    objectStr = "HA = parallelHybridAutomaton(comp,iBinds);";
    
    % likely redundant:
%     if all(strcmp(cellfun(@(comp) comp.inputs.name,components),'uDummy'))
%         % if no inputs in system, CORA needs an input
%         % generate a global input on the first component
%         % this will have no effect as the input matrices should all be 0
%         componentStr = componentStr + ...
%             "% no inputs given, global input for CORA, no effect as B = 0" + ...
%             newline + "iBinds{1} = [0,1];" + aux_newlines(2);
%     end
end

% concatenate all strings
HA = automatonStr + componentStr + objectStr + aux_newlines(2) + "end";

end


% Auxiliary functions -----------------------------------------------------

% main auxiliary functions
function infoStr = aux_infoStr(components)
    
% number of components
nrComp = length(components);
% flat or parallel conversion
isFlatHA = nrComp == 1;

% init string
infoStr = "";

% loop over all components
for iComp = 1:nrComp

    % read out i-th component for quicker access
    comp = components{iComp};
    % read out name of component
    compInfoStr = "% Component " + int2str(iComp) + " (" + comp.name + "):" + newline;
    
    % gather state names in string array
    if ~isempty(comp.states)
        stateStr = "%  state x := [" + strjoin([comp.states.name],"; ") + "]" + newline;
    else
        stateStr = "%  state x := [];" + newline;
    end
    
    % gather input names in string array
    if ~isempty(comp.inputs)
        inputStr = "%  input u := [" + strjoin([comp.inputs.name],"; ") + "]" + newline;
    else
        stateStr = "%  input u := [];" + newline;
    end

    % gather output names in string array
    if isFlatHA
        % flat: accesses outputsGlobal
        if ~isempty(comp.outputsGlobal)
            outputStr = "%  output y := [" + strjoin([comp.outputsGlobal.name],"; ") + "]" + aux_newlines(2);
        else
            outputStr = newline;
        end
    else
        % parallel: accesses outputsLocal
        if ~isempty(comp.outputsLocal)
            outputStr = "%  output y := [" + strjoin([comp.outputsLocal.name],"; ") + "]" + aux_newlines(2);
        else
            outputStr = newline;
        end
    end
    
    infoStr = infoStr + compInfoStr + stateStr + inputStr + outputStr;
end

end

function dynamicsStr = aux_dynamicsStr(loc,iLoc,comp,iComp,isFlatHA,resultpath,functionName)

% write original equation as a comment
dynamicsCommFlow = aux_text2comment("flow equation:" + newline + loc.Flow.Text) + newline;
dynamicsCommOutput = aux_text2comment("output equation:" + newline + loc.Invariant.Text_output) + newline;

% linear of nonlinear flow/output?
if isfield(loc.Flow,'A')
    % only enters here if both flow equation and output equation are linear

    % read out state equation of linear system
    linSys_A = printMatrixConverter(loc.Flow.A);
    linSys_AStr = "dynA = ..." + newline + linSys_A + ";" + newline;
    linSys_B = printMatrixConverter(loc.Flow.B);
    linSys_BStr = "dynB = ..." + newline + linSys_B + ";" + newline;
    linSys_c = printMatrixConverter(loc.Flow.c);
    linSys_cStr = "dync = ..." + newline + linSys_c + ";" + newline;
    
    % outputs (which can only be inputs to other components) are
    % defined in the invariant of the location
    if (isFlatHA && isempty(comp.outputsGlobal)) || ...
        (~isFlatHA && isempty(comp.outputsLocal) && isempty(comp.outputsGlobal))
        % no outputs -> only state equation x' = Ax + Bu + c
        dynamicsStr = dynamicsCommFlow + linSys_AStr + linSys_BStr + linSys_cStr + ...
            "dynamics = linearSys(dynA, dynB, dync);" + aux_newlines(2);
    else
        % read out output equation of linear system
        linSys_C = printMatrixConverter(loc.Flow.C);
        linSys_CStr = "dynC = ..." + newline + linSys_C + ";" + newline;
        linSys_D = printMatrixConverter(loc.Flow.D);
        linSys_DStr = "dynD = ..." + newline + linSys_D + ";" + newline;
        linSys_k = printMatrixConverter(loc.Flow.k);
        linSys_kStr = "dynk = ..." + newline + linSys_k + ";" + newline;
        % two equations: x' = Ax + Bu + c and y = Cx + Du + k
        dynamicsStr = dynamicsCommFlow + linSys_AStr + linSys_BStr + linSys_cStr + ...
            newline + dynamicsCommOutput + linSys_CStr + linSys_DStr + linSys_kStr + ...
            "dynamics = linearSys(dynA, dynB, dync, dynC, dynD, dynk);" + aux_newlines(2);
    end
else
    % choose name for dynamics function
    if isFlatHA
        % simpler name for flat automata
        nonlinName = sprintf("%s_Loc%d_FlowEq",functionName,iLoc);
        nonlinName_out = sprintf("%s_Loc%d_OutputEq",functionName,iLoc);
    else
        nonlinName = sprintf("%s_Comp%d_Loc%d_FlowEq",functionName,iComp,iLoc);
        nonlinName_out = sprintf("%s_Comp%d_Loc%d_OutputEq",functionName,iComp,iLoc);
    end
    
    % number of states x and inputs u of the flow and output equation
    statedims = num2str(length(comp.states));
    inputdims = num2str(length(comp.inputs));
    
    % write file for nonlinear flow
    printDynamicsFile(resultpath,nonlinName,loc.Flow.FormalEqs,"flow");
    
    if isFlatHA || strlength(loc.Invariant.Text_output) == 0
        % no output equation (always in conversion to flat HA)
        dynamicsStr = dynamicsCommFlow + "dynamics = nonlinearSys(@" + ...
            nonlinName + "," + statedims + "," + inputdims + "); " + newline;

    else
        % output equation -> write file
        printDynamicsFile(resultpath,nonlinName_out,loc.Flow.FormalEqs_out,"flow");

        % number of outputs
        outputdims = num2str(length(comp.outputsLocal));
        
        % write string for nonlinear flow equation
        dynamicsStr = dynamicsCommFlow + dynamicsCommOutput + ...
            "dynamics = nonlinearSys(@" + nonlinName + "," + statedims + ...
            "," + inputdims + ",..." + newline + "    " + ...
            "@" + nonlinName_out + "," + outputdims + "); " + aux_newlines(2);
    end
end

end

function invariantStr = aux_invariantStr(loc,allTrans,isFlatHA)

% Get information for Invariant
InvText = loc.Invariant.Text;
if isa(loc.Invariant.set,'fullspace')
    n = length(loc.Flow.expressions);
    [str1,str2] = aux_fullspaceString(n,allTrans,isFlatHA);
elseif isa(loc.Invariant.set,'polytope')
    [str1,str2] = aux_polytopeString(loc.Invariant.set);
elseif isa(loc.Invariant.set,'levelSet')
    [str1,str2] = aux_levelSetString(loc.Invariant.set);
else
    throw(CORAerror('CORA:converterIssue',...
        'Invariant has to be either empty, a polytope, or a levelSet'));
end

% Write String for Invariant
invariantComm = aux_text2comment("invariant equation:" + newline + InvText) + newline;
invariantStr = invariantComm + str1 + "inv = " + str2 + aux_newlines(2);

end

function resetStr = aux_resetStr(trans,iComp,iLoc,iTrans,isFlatHA,resultpath,functionName)

% linear or nonlinear reset function?
if isfield(trans.reset,'A')
    % linear reset
    reset_A = printMatrixConverter(trans.reset.A);
    reset_AStr = "resetA = ..." + newline + reset_A + ";" + newline;
    if isfield(trans.reset,'B')
        reset_B = printMatrixConverter(trans.reset.B);
        reset_BStr = "resetB = ..." + newline + reset_B + ";" + newline;
    end
    reset_c = printMatrixConverter(trans.reset.c);
    reset_cStr = "resetc = ..." + newline + reset_c + ";" + newline;
    
    % Write Reset String
    tranResetText = trans.reset.Text;
    if tranResetText == ""
        tranResetText = "no reset equation given";
    end
    resetComm = aux_text2comment("reset equation:" + newline + tranResetText) + newline;
    if isfield(trans.reset,'B')
        resetStr = resetComm + reset_AStr + reset_BStr + reset_cStr + ...
            "reset = struct('A', resetA, 'B', resetB, 'c', resetc);" + aux_newlines(2);
    else
        resetStr = resetComm + reset_AStr + reset_cStr + ...
            "reset = struct('A', resetA, 'c', resetc);" + aux_newlines(2);
    end
else
    % nonlinear reset
                   
    % choose name of function
    if isFlatHA
        % flat hybrid automaton
        resetFuncName = sprintf("%s_Loc%d_Trans%d_ResetEq%d",...
            functionName,iLoc,iTrans);
    else
        % parallel hybrid automaton
        resetFuncName = sprintf("%s_Comp%d_Loc%d_Trans%d_ResetEq%d",...
            functionName,iComp,iLoc,iTrans);
    end

    % generate function file for nonlinear reset
    printDynamicsFile(resultpath,resetFuncName,trans.reset.FormalEqs,"reset");
    tranResetText = trans.reset.Text;
    if tranResetText == ""
        tranResetText = "no reset equation given";
    end
    
    % Write Reset String
    resetComm = aux_text2comment("reset equation:" + ...
        newline + tranResetText) + newline;
    resetStr = resetComm + "reset = struct('f', @" + ...
            resetFuncName + ");" + aux_newlines(2);
end

end

function guardStr = aux_guardStr(trans,loc)

% general information about guard set
tranGuardText = trans.guard.Text;
if tranGuardText == ""
    tranGuardText = "no guard set given -> instant transition";
end

% different set representations for guard set
if isa(trans.guard.set,'fullspace')
    % fullspace -> immediate transition (no text in spaceex)
    str1 = "";
    % init invariant using state dimension of dynamics in given location
    str2 = "fullspace(" + length(loc.Flow.expressions) + ");";

elseif isa(trans.guard.set,'polytope')
    % intersect guard with invariant
    try
        G = loc.Invariant.set & trans.guard.set;                    
    catch
        G = trans.guard.set;
    end

    % check if guard can be represented as hyperplane
    res = false;
    if isa(G,'polytope')
        G = compact_(G,'all',1e-9);
        [res,hyp] = representsa_(G,'conHyperplane',eps);
    end
    
    if isa(G,'levelSet')
        [str1,str2] = aux_levelSetString(G);
    elseif isa(G,'polytope')
        if res
            [str1,str2] = aux_conHyperplaneString(hyp);
        else
            [str1,str2] = aux_polytopeString(G);
        end
    else
        throw(CORAerror('CORA:converterIssue',...
            'Unexpected set representation in conversion of guard set.'));
    end

elseif isa(trans.guard.set,'levelSet')
    % intersect guard with invariant
    try
        G = loc.Invariant.set & trans.guard.set;                    
    catch
        G = trans.guard.set;
    end
    [str1,str2] = aux_levelSetString(G);

else
    % throw error in case unexpected set representation comes up
    throw(CORAerror('CORA:converterIssue',...
        ['Guard set has to be empty, a conHyperplane, '...
        'a polytope or a levelSet.'])); 
end

% write guard string
guardComm = aux_text2comment("guard equation:" + newline + tranGuardText) + newline;
guardStr = guardComm + str1 + "guard = " + str2 + aux_newlines(2);

end

function inputBindsStr = aux_inputBindsStr(comp,iComp,components)

% general information about input binds
parallelComm = aux_text2comment("composition: hybrid automaton and input binds") + ...
    newline + "comp(" + num2str(iComp) + ") = hybridAutomaton(loc);" + aux_newlines(2);

% syntax of inputBinds: mx2 array, where
%   m - number of input arguments to current component
% with (m,1) being either 0 (global) or the component of origin
% and (m,2) being the number of the global input (if (m,1) = 0)
%           or the number of output of the component of origin

if strcmp(comp.inputs(1).name,'uDummy')
    % only dummy input generated by conversion... use first global
    % input as input bind (no effect as all input matrices are 0)
    inputBindsStr = "% only dummy input" + newline + ...
        "iBinds{" + iComp + "} = [0 1];" + aux_newlines(2);

else
    % actual (meaningful) inputs given

    inputBindsStr = "% input names: " + strjoin([comp.inputs.name],", ") ...
        + newline + "iBinds{" + iComp + "} = [[";

    % loop over all inputs
    for inp=1:length(comp.inputs)
        if inp > 1
            inputBindsStr = inputBindsStr + ";[";
        end
        inputName = comp.inputs(inp).name;
        % default: input comes from composed system ('global')
        origin = 0;
        outputOfOrigin = 1;
        % search for inputName in all other components
        for c=1:length(components)
            idx = strcmp([components{c}.outputsLocal.name],inputName);
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
        inputBindsStr = inputBindsStr + num2str(origin) + "," +...
            num2str(outputOfOrigin) + "]";
    end
    inputBindsStr = inputBindsStr + "];" + aux_newlines(2);
end

inputBindsStr = parallelComm + inputBindsStr;

end


% writing sets as string
function [str1,str2] = aux_fullspaceString(n,allTrans,isFlatHA)
% in principle, we interpret an undefined invariant in SpaceEx (no text) as
% the entire n-dim. space, corresponding to a fullspace object in CORA; the
% reachability analysis requires to leave the invariant in order to proceed
% to the transition to another location; therefore, we need to 'cut' the
% invariant according to the outgoing transitions of the location

if isempty(allTrans) || n == 0 || (~isFlatHA && ...
        all(arrayfun(@(x) ~strcmp(x,''),[allTrans.label],'UniformOutput',true)) )
    % conditions for fullspace invariant:
    % 1. no outgoing transitions -> fullspace as invariant
    % 2. no flow/invariant given -> dummy fullspace (since n=0)
    % 3. all outgoing transitions have a synchronization label -> this
    %    component has to wait for other components to leave their
    %    invariant, which will trigger the transition
    % additionally, if 
    str1 = "";
    str2 = "fullspace(" + n + ");";
else
    % instantiate fullspace object
    inv = fullspace(n);
    % loop over all outgoing transitions and shrink the invariant by the
    % complements of each guard set
    for t=1:length(allTrans)
        % skip guards with synchronization label (will be computed
        % on-the-fly, if the entire transition is active)
        if isFlatHA || strcmp(allTrans(t).label,'')
            % read out guard set
            guard = allTrans(t).guard.set;
            % shift by a small bit if it's a polytope with a strict
            % inequality -> necessary to result in empty invariants for
            % pairs of guard sets ax <= b & ax > b
            if isa(guard,'polytope') && ...
                  ( contains(allTrans(t).guard.Text,"<") || ...
                    contains(allTrans(t).guard.Text,">") ) && ...
                 ~( contains(allTrans(t).guard.Text,">=") || ...
                    contains(allTrans(t).guard.Text,"<=") )
                guard = guard + 1e-12*ones(n,1);
            end
    
            % try to compute complement of guard set
            if isa(guard,'polytope') || isa(guard,'levelSet') ...
                    || isa(guard,'fullspace')
                temp = ~guard;
            else
                % complement cannot be computed
                throw(CORAerror('CORA:converterIssue',...
                    'Issue in conversion of empty invariant and non-empty guard sets.'));
            end
    
            % intersect with invariant
            inv = and_(inv,temp,'exact');
    
            % break if invariant has become empty
            if isa(inv,'emptySet')
                break
            end
        end
    end
    % write to string (either polytope or levelSet)
    if isa(inv,'emptySet') || (isa(inv,'polytope') && isempty(inv))
        % resulting invariant is the empty set -> the location will be
        % exited even before the first step of the reachability analysis
        str1 = ""; str2 = "emptySet(" + n + ");";
    elseif isa(inv,'polytope')
        [str1,str2] = aux_polytopeString(inv);
    elseif isa(inv,'levelSet')
        [str1,str2] = aux_levelSetString(inv);
    else
        throw(CORAerror('CORA:converterIssue',...
            'Issue in conversion of empty invariant.'));
    end
end

end

function [str1,str2] = aux_polytopeString(set)
% generates the string that constructs the polytope

    A = printMatrixConverter(set.A);
    AStr = "P_A = ..." + newline + A + ";" + newline;
    b = printMatrixConverter(set.b);
    bStr = "P_b = ..." + newline + b + ";" + newline;
    if ~isempty(set.Ae)
        % if invariant is a polyhedron, add additional parameters
        Ae = printMatrixConverter(set.Ae);
        AeStr = "P_Ae = ..." + newline + Ae + ";" + newline;
        be = printMatrixConverter(set.be);
        beStr = "P_be = ..." + newline + be + ";" + newline;
        str2 = 'polytope(P_A,P_b,P_Ae,P_be);';
    else
        AeStr = "";
        beStr = "";
        str2 = 'polytope(P_A,P_b);';
    end
    str1 = AStr + bStr + AeStr + beStr + newline;
end

function [str1,str2] = aux_levelSetString(set)
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
    str1 = varStr + newline + eqStr + newline + compStr + aux_newlines(2);
    str2 = "levelSet(eq,vars,compOp);";

end

function [str1,str2] = aux_conHyperplaneString(set)
% generates the string that constructs the conHyperplane object

    % hyperplane equation c*x = d
    c = printMatrixConverter(set.a');
    cStr = "c = " + c + ";" + newline;
    dStr = "d = " + num2str(set.b) + ";";
    
    % inequality constraints C*x <= D
    if ~isempty(set.C)
        C = printMatrixConverter(set.C);
        Cstr = "C = ..." + newline + C + ";" + newline;
        D = printMatrixConverter(set.d);
        Dstr = "D = " + D + ";";
        
        str1 = cStr + dStr + Cstr + Dstr + aux_newlines(2);
        str2 = "conHyperplane(c,d,C,D);";
    else
        str1 = cStr + dStr + aux_newlines(2);
        str2 = "conHyperplane(c,d);";
    end
end


% helper functions for formatting
function str = aux_padComment(comment,maxLineLength)
%pads comment left & right with dashes to desired length and prefixes "%"

    if nargin < 2
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
