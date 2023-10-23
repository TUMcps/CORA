function [isLinear,A,B,c,formalEqs] = ...
    eq2linSysFlat(exprNames,exprs,states,inputs,outputsLocal,map,parseMode)
% eq2linSysFlat - Converts string equations to a LTI system specification
%    (if possible). Returns formatted nonlinear equations otherwise.
%    Requires names of state and input variables as parameter.
%    Requires names of constants and values to substitute.
%    WARNING: could produce incorrect results, if any variables are named 
%       "xL<number>R" or "uL<number>R
%
% Syntax:
%    [isLinear,A,B,c,formalEqs] = ...
%       eq2linSysFlat(exprNames,exprs,states,inputs,outputsLocal,map,parseMode)
%
% Inputs:
%    exprNames (symbolic) - Ordered list specifing which expression
%                           defines which state, used to correctly order
%                           the equations needed for a nonlinear system
%    exprs (string) - flow equation in SX form
%              rough format := <eq>("&" <eq>)*
%                      <eq> := <state>' "==" <expr(states,inputs,constants)>
%              (permitting '=', '==', or ':=' for assignment)
%    states (struct) - states of LTI system
%    inputs (struct) - input signals
%    outputsLocal - ???
%    map - ???
%    parseMode - affects default outputs ('flow'(default), 'assignment')
%
% Outputs:
%    isLinear - equations successfully linearized? (logical)
%    A - state matrix of LTI system
%    B - input matrix of LTI system
%    c - constant offset of LTI system
%    formalEqs - equations formatted to paste into a matlab dynamics file
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       ???
% Written:       ---
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin < 7
    parseMode = 'flow';
elseif ~strcmp(parseMode,'flow') && ~strcmp(parseMode,'assignment')
    warning('Unknown mode "%s", using default:"flow"',parseMode);
    parseMode = 'flow';
end

%create symbolic variables for states
%(with L/R brackets to later replace with regular brackets)
numStates = length(states);

stateNames = strings(numStates,1);
varnames = cell(numStates,1);
for i=1:numStates
     stateNames(i) = states(i).name;
     varnames{i} = strcat('xL',num2str(i),'R');
end
x = sym(varnames);

%create symbolic variables for inputs
numInputs = length(inputs);

inputNames = strings(numInputs,1);
varnames = cell(numInputs,1);
for i=1:numInputs
    inputNames(i) = inputs(i).name;
    varnames{i} = strcat('uL',num2str(i),'R');
end
u = sym(varnames);

%create symbolic variables for local outputs
numOutputs = length(outputsLocal);

outputNames = strings(numOutputs,1);
varnames = cell(numOutputs,1);
for i=1:numOutputs
    outputNames(i) = outputsLocal(i).name;
    varnames{i} = strcat('yL',num2str(i),'R');
end
y = sym(varnames);

% substitute states into equation
exprs = applySymMapping(exprs,stateNames,x);

% substitute inputs into equation
exprs = applySymMapping(exprs,inputNames,u);

% substitute outputsLocal into equation
if ~isempty(y)
    j = 1;
    for i=1:length(stateNames)
        if strcmp(map,stateNames(i))
            y_x(j) = x(i);
            j = j+1;
        end
    end
    exprs = applySymMapping(exprs,outputNames,y_x);
end

%check whether variables except x & u & y remain
isSubset = all(ismember(symvar(exprs), [x;u;y]));
if ~isSubset
    %print detailed error message
    diff = setdiff(symvar(exprs), [x;u;y]);
    varstr = char(diff(1));
    if(length(diff) > 1)
        varstr = [varstr sprintf(', %s',diff(2:end))];
    end
    throw(CORAerror('CORA:converterIssue',...
        ['Error while parsing ' char(parseMode) ': unknown variables [' varstr ']']));
end

% assign exprs to symbolic function, such that f(i) = xi'
f = sym('f',[numStates,1]);
for i = 1:numStates
    eqIdx = 0;
    for j = 1:length(exprs)
        if strcmp(states(i).name, exprNames{j})
            eqIdx = j;
            break;
        end
    end
    if eqIdx ~= 0
        f(i,1) = exprs(j);
    else
        if strcmp(parseMode,'assignment')
            f(i,1) = x(i);
            % ommitting assignments allowed per SX specification
        else % assuming default mode 'flow'
            f(i,1) = sym(0);
            warning("No flow specified for state %s!\n" +...
            "assuming %s' = 0...",states(i).name, states(i).name);
        end
    end
end

%computing linear dependencies
A_sym = jacobian(f,x);
B_sym = jacobian(f,u);

%computing constant components, by setting x = u = 0
f_const = subs(f, [x;u], zeros(numStates + numInputs,1));
%evaluate as column vector
c_sym = f_const(1:numStates);

try
    A = double(A_sym);
    B = double(B_sym);
    c = double(c_sym);
    isLinear = true;
catch
    %Could not linerarize system.
    A = [];
    B = [];
    c = [];
    isLinear = false;

    % Dear future person, who is tasked to implement linear systems with
    % variable parameters:  
    %   When A is allowed to depend on parameters, simply check whether
    %   symvar(A_sym) contains only parameter variables. Then A_sym(p)
    %   could be rewritten as a matlab function, or something like that.
end

% Format equations for use in matlab function
% string operations may fail in versions older than R2016b
eqs = strings(1,numStates);
for i = 1:numStates
    eq = string(f(i));
    eq = strrep(eq,'L','(');
    eq = strrep(eq,'R',')');
    eqs(i) = sprintf("dx(%d,1) = %s;\n",i,eq);
end

formalEqs = join(eqs,newline);

% ------------------------------ END OF CODE ------------------------------
