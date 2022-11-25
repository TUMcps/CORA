function [C, D, k] = outputMatrices(EqExprs, states, inputs, outputs)
% Converts string equalities defining the output of a location
% to the C, D and k matrices required for linearSys instantiation
% WARNING: could produce incorrect results, if any variables are named 
% "xL<number>R", "uL<number>R" or "yL<number>R"
%
% Remark: broadly copied syntax from eq2polytope.m
%
%
% Syntax:  
%    [C, D, k] = outputMatrices(EqExprs, listOfVar, states, inputs, outputs)
%
% Inputs:
%    EqExprs   - equations defining outputs
%      equation: string - rough format: <eq> ("&" <eq>)*
%                  <eq> = <expr(states,constants)> <op> <expr(states,constants)>
%                  <expr>: linear expressions only
%                  <op> = "<"|">"|"<="|">="|==
%    states    - all state variables
%    inputs    - all input variables
%    outputs   - all output variables
%
% Outputs:
%    C - output matrix
%    D - throughput matrix
%    k - output offset
%
% Other m-files required: none
% Subfunctions: ADD IF NECESSARY
% MAT-files required: none
%
% See also: none

% Author: Mark Wetzlinger
% Written: 26-December-2018 
% Last update: --- 
% Last revision: ---

%------------- BEGIN CODE --------------

% if no outputs, then no equation y = Cx + Du + k
%   hence, return empty matrices (different handling in data2parallelHA.m)
if isempty(outputs) || isempty(EqExprs)
    C = [];
    D = [];
    k = [];
    return
end

% size init
sizeC = length(states);
sizeD = length(inputs);
sizek = length(outputs);
C = zeros(sizek,sizeC);
D = zeros(sizek,sizeD);
k = zeros(sizek,1);

% create symbolic variables for all variables
numvars = length(states) + length(inputs) + length(outputs);
% remark: cannot take listOfVar because of uDummy
varnames = cell(numvars,1);
listOfVarNames = strings(numvars,1);

for i=1:sizeC
    listOfVarNames(i) = states(i).name;
    varnames{i} = strcat('xL',num2str(i),'R');
end

j = i+1;
for i=1:sizeD
    listOfVarNames(j) = inputs(i).name;
    varnames{j} = strcat('uL',num2str(i),'R');
    j = j+1;
end

for i=1:sizek
    listOfVarNames(j) = outputs(i).name;
    varnames{j} = strcat('yL',num2str(i),'R');
    j = j+1;
end

x = sym(varnames);

% substitute variables into equations
EqExprs = applySymMapping(EqExprs,listOfVarNames,x);

% check whether any expressions include non-state variables
% (compute a logical index for this)
eqValidIdx = false(size(EqExprs));
for i=1:length(eqValidIdx)
    eqValidIdx(i) = all(ismember(symvar(EqExprs(i)), x));
end

% if any expressions fail the test, print warnings then remove them
if(~all(eqValidIdx))
    % print warnings
    for i=1:length(eqValidIdx)
        if ~eqValidIdx(i)
            % detect, which variables caused the error
            diff = setdiff(symvar(EqExprs(i)), x);
            % print detailed warning message
            varstr = "(" + string(diff(1));
            if(length(diff) > 1)
                varstr = varstr + sprintf(", %s",diff(2:end));
            end
            varstr = varstr + ")";
            warning("A CONDITION CONTAINS NON-STATE VARIABLES %s AND IS IGNORED!"+...
                newline + "Condition: ""%s == 0"" (after arithmetic transformation)",...
                varstr,string(EqExprs(i)));
        end
    end
    % remove invalid expressions from arrays
    EqExprs = EqExprs(eqValidIdx);
end

% computing linear dependencies
C_sym = jacobian(EqExprs,x);

% computing constant components, by setting x = 0
k = subs(EqExprs, x, zeros(length(x),1));
k = double(k) * -1;

% compute C matrix
% assumption: always at least one state given -> no length check
% find correct sequence of to-be-written rows of C
% for the case when: y2 = x1 and y1 = x2 (rowsequence -> [2 1])
% equations get switched (only if more than one output)
for r=1:size(C_sym,1)
    temp = find(C_sym(r,:));
    posinC = temp(find(C_sym(r,:)) > length(states) + length(inputs));
    rowsequence(r) = posinC - length(states) - length(inputs);
end

xIdx = strfind(varnames,'x');
for row=1:size(C_sym,1)
    CIdx = 1;
    for i=1:numvars
        if xIdx{i}
            C(rowsequence(row),CIdx) = - C_sym(row,i);
            CIdx = CIdx + 1;
        end
    end
end


% compute D matrix
if ~isempty(inputs) && ~strcmp(inputs(1).name,'uDummy')
    % D ~= 0
    % find all entries in varnames which contain a 'u' -> input!
    uIdx = strfind(varnames,'u');
    for row=1:size(C_sym,1)
        % counter for position in D matrix
        DIdx = 1;
        for i=1:numvars
            if uIdx{i}
                % index in varnames corresponds to input and to position in D
                D(rowsequence(row),DIdx) = -C_sym(row,i);
                DIdx = DIdx + 1;
            end
        end
    end
else
    D = 0;
end

end

%------------- END OF CODE -------------