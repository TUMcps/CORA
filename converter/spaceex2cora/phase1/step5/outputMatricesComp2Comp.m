function [isLin_out,C,D,k,eqs_out] = outputMatricesComp2Comp(states,inputs,outputs,inv)
% outputMatricesComp2Comp - instantiates output matrices C, D, and output
%    vector k for output matrices of a component in a parallel hybrid
%    automaton
%
% Syntax:
%    [isLin_out,C,D,k,eqs_out] = outputMatricesComp2Comp(states,inputs,outputs,inv)
%
% Inputs:
%    states - all state variables
%    inputs - all input variables
%    outputs - all output variables
%    inv - invariant (containing output equations)
%
% Outputs:
%    isLin_out - flag whether all output equations are linear
%    C - output matrix
%    D - feedthrough matrix
%    k - output offset
%    eqs_out - string containing text for nonlinear output equations
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       25-June-2022
% Last update:   ---
% Last revision: 12-January-2023 (MW, major fixes, nonlinear output eqs)

% ------------------------------ BEGIN CODE -------------------------------

% if no outputs, then no equation y = Cx + Du + k
% hence, return empty matrices (different handling in data2parallelHA.m)
if isempty(outputs)
    C = []; D = []; k = []; isLin_out = true; eqs_out = "";
    return
end

% size init
sizeC = length(states);
sizeD = length(inputs);
sizek = length(outputs);
C = zeros(sizek,sizeC);
D = zeros(sizek,sizeD);
k = zeros(sizek,1);

% identity matrices
Ix = eye(length(states));
Iu = eye(length(inputs));

% nonlinear output equations
eqs_out = "";

% convert to matrices C, D, and k until nonlinear equation found
isLin_out = true;

% map all states into one state vector for check
sym_state = sym('x',[length(states),1]);
sym_input = sym('u',[length(inputs),1]);

% loop over all outputs
for i=1:length(outputs)
    % outputs can be defined with/without an explicit output equation:
    % - implicitly: only linked to other base component in SpaceEx model
    % - explicitly: output equation in invariant of location

    % see if there is an output equation for the i-th output in the
    % invariant; left side should be a vector of the output variable names
    idxOutputEq = logical(inv.exprLeft - outputs(i).name == 0);

    % no explicit output equation -> construct an identity output
    if ~any(idxOutputEq)
        % see to which state the i-th output corresponds
        outputIdx = [states.name] == outputs(i).name;
        % add entry in C matrix
        C(i,outputIdx) = 1;
        % go to next output
        continue
    end


    % replace chosen variables in output equation
    outputEq = inv.exprRight(i);
    for j=1:length(states)
        outputEq = subs(outputEq,states(j).name,sym_state(j));
    end
    for j=1:length(inputs)
        outputEq = subs(outputEq,inputs(j).name,sym_input(j));
    end
    % convert to string
    outputEqStr = string(outputEq);
    % add brackets
    for j=1:length(states)
        outputEqStr = strrep(outputEqStr,"x" + j,"x(" + j + ")");
    end
    for j=1:length(inputs)
        outputEqStr = strrep(outputEqStr,"u" + j,"u(" + j + ")");
    end

    % list in output equations: we use dx to comply with the standard
    % rewriting into functions known from differential equations
    eqs_out = eqs_out + "dx(" + i + ",1) = " + outputEqStr + ";" + newline;

    % instantiate function handle (for isFuncLinear call)
    fun = eval("@(x,u) " + outputEqStr + ";");
    if ~isLin_out || ~isFuncLinear(fun)
        isLin_out = false; continue
    end

    % output matrix C
    temp = subs(outputEq,sym_input,zeros(length(inputs),1));
    for j=1:length(states)
        % substitute 1 for j-th state and 0 for all other states to
        % determine coefficient in C matrix 
        C(i,j) = double(subs(temp,sym_state,Ix(:,j)));
    end

    temp = subs(outputEq,sym_state,zeros(length(states),1));
    % feedthrough matrix D
    for j=1:length(inputs)
        % substitute 1 for j-th state and 0 for all other states to
        % determine coefficient in C matrix
        D(i,j) = double(subs(temp,sym_input,Iu(:,j)));
    end

    % constant offset k: all that remains
    temp = subs(outputEq,sym_state,zeros(length(states),1));
    k(i) = double(subs(temp,sym_input,zeros(length(inputs),1)));

end

% ------------------------------ END OF CODE ------------------------------
