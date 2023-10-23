function [isLin_out,C,D,k,eqs_out] = ...
    outputMatrices(EqExprs,states,inputs,outputsLocal,map,outputsGlobal)
% outputMatrices - Converts string equalities defining the output of a
%    location to the C, D matrices and k vector required for linearSys
%    WARNING: could produce incorrect results, if any variables are named 
%             "xL<number>R", "uL<number>R", or "yL<number>R"
%
% Remark: broadly copied syntax from eq2polytope.m (now eq2set.m)
%
% Syntax:
%    [isLin_out,C,D,k,eqs_out] = ...
%       outputMatrices(EqExprs,states,inputs,outputsLocal,map,outputsGlobal)
%
% Inputs:
%    EqExprs - equations defining outputs, rough format: <eq> ("&" <eq>)*
%              <eq> = <expr(states,constants)> <op> <expr(states,constants)>
%              <expr>: linear expressions only
%              <op> = "<"|">"|"<="|">="|==
%    states - all state variables
%    inputs - all input variables
%    outputsLocal - output variables that can be resolved to states/inputs
%    map - where outputsLocal maps to
%    outputsGlobal - global outputs (for which the output equations are
%                    computed)
%
% Outputs:
%    isLin_out - true/false whether output equation linear
%    C - output matrix
%    D - feedthrough matrix
%    k - output offset
%    eqs_out - output equation in text form
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       26-December-2018 
% Last update:   ---
% Last revision: 22-January-2023 (MW, rewrite & additional output arguments)

% ------------------------------ BEGIN CODE -------------------------------

% if no outputs, then no equation y = Cx + Du + k
% hence, return empty matrices (see processing in data2parallelHA.m)
if isempty(outputsGlobal) || isempty(EqExprs)
    C = []; D = []; k = []; isLin_out = true; eqs_out = "";
    return
end

% momentarily not done...
eqs_out = "";

% number of states, inputs, and global outputs
nrStates = length(states);
nrInputs = length(inputs);
nrOutputs = length(outputsGlobal);
% initialize matrices C, D, and k
C = zeros(nrOutputs,nrStates);
D = zeros(nrOutputs,nrInputs);
k = zeros(nrOutputs,1);

% create symbolic variables for states
stateNames = strings(nrStates,1);
varnames = cell(nrStates,1);
for i=1:nrStates
     stateNames(i) = states(i).name;
     varnames{i} = strcat('xL',num2str(i),'R');
end
x = sym(varnames);

% create symbolic variables for inputs
inputNames = strings(nrInputs,1);
varnames = cell(nrInputs,1);
for i=1:nrInputs
    inputNames(i) = inputs(i).name;
    varnames{i} = strcat('uL',num2str(i),'R');
end
u = sym(varnames);

% create symbolic variables for global outputs
outputNames = strings(nrOutputs,1);
varnames = cell(nrOutputs,1);
for i=1:nrOutputs
    outputNames(i) = outputsGlobal(i).name;
    varnames{i} = strcat('yL',num2str(i),'R');
end
yGlobal = sym(varnames);

% all variables
allVars = [x;u;yGlobal];

% copy EqExprs to keep original
EqExprsSym = EqExprs;

% substitute local outputs into equations
if ~isempty(outputsLocal)
    yLocal = sym([]);
    % map to states?
    for i=1:length(states)
        if strcmp(map,stateNames(i))
            yLocal(end+1,1) = x(i);
        end
    end
    % map to inputs?
    for i=1:length(inputs)
        if strcmp(map,inputNames(i))
            yLocal(end+1,1) = u(i);
        end
    end
    EqExprsSym = applySymMapping(EqExprsSym,[outputsLocal.name],yLocal);
end

% substitute remaining variables into equations
EqExprsSym = applySymMapping(EqExprsSym,...
    [stateNames;inputNames;outputNames],allVars);

% logical index whether any expressions include non-state/input variables
eqValidIdx = false(size(EqExprsSym));
for i=1:length(eqValidIdx)
    eqValidIdx(i) = all(ismember(symvar(EqExprsSym(i)), allVars));
end

% if any expression fails the test, print warnings and remove them
if(~all(eqValidIdx))
    % print warnings
    for i=1:length(eqValidIdx)
        if ~eqValidIdx(i)
            % detect, which variables caused the error
            diff = setdiff(symvar(EqExprsSym(i)), x);
            % print detailed warning message
            varstr = "(" + string(diff(1));
            if(length(diff) > 1)
                varstr = varstr + sprintf(", %s",diff(2:end));
            end
            varstr = varstr + ")";
            warning("A CONDITION CONTAINS NON-STATE/INPUT VARIABLES %s AND IS IGNORED!"+...
                newline + "Condition: ""%s == 0"" (after arithmetic transformation)",...
                varstr,string(EqExprsSym(i)));
        end
    end
    % remove invalid expressions from arrays
    EqExprsSym = EqExprsSym(eqValidIdx);
end

% remove all equations that do not contain any global output variable
isOutputEq = false(size(EqExprsSym));
for i=1:length(isOutputEq)
    for j=1:length(outputNames)
        if has(EqExprsSym(i),yGlobal(j))
            isOutputEq(i) = true; break
        end
    end
end
% remove all non-output equation
EqExprsSym = EqExprsSym(isOutputEq);

% computing constant components, by setting x = 0
k = subs(EqExprsSym, allVars, zeros(length(allVars),1));
k = double(k) * -1;

% compute C and D matrices via Jacobian
C_sym = jacobian(EqExprsSym,x);
% compute D matrix
D_sym = jacobian(EqExprsSym,u);


% convert Jacobian to double arrays
try
    % switch sign since equations are given in the form, e.g.,
    %    y1 - x1 (= 0)  <=>  y1 = x1
    C = -double(C_sym);
    D = -double(D_sym);
    % text for comment in model file
    EqStr = cell(0);
    for i=1:length(EqExprs)
        if isOutputEq(i)
            % substitute 0 for global outputs and invert sign
            EqStr{end+1,1} = char("y" + i + " = " ...
                + string(-subs(EqExprs(i),outputNames,zeros(nrOutputs,1))));
        end
    end
    % concatenate text of individual output equations
    eqs_out = strjoin(EqStr," & ");
    % output equations are linear
    isLin_out = true;
catch
    isLin_out = false; C = []; D = []; k = [];
end

% ------------------------------ END OF CODE ------------------------------
