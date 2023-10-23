function [isLinear,A,B,c,formalEqs,containsInput] = ...
    eq2linSys(exprNames,exprs,states,inputs,varargin)
% eq2linSys - Converts string equations to a LTI system specification
%    (if possible). Returns formatted nonlinear equations otherwise.
%    Requires names of state and input variables as parameter.
%    Requires names of constants and values to substitute.
%    WARNING: could produce incorrect results, if any variables are named 
%       "xL<number>R" or "uL<number>R
%
% Syntax:
%    [isLinear,A,B,c,formalEqs,containsInput] = ...
%       eq2linSys(exprNames,exprs,states,inputs,parseMode)
%
% Inputs:
%    exprNames (symbolic) - Ordered list specifing which expression
%                           defines which state, used to correctly order
%                           the equations needed for a non-linear system
%    exprs (string) - flow equation in SX form
%              rough format := <eq>("&" <eq>)*
%                      <eq> := <state>' "==" <expr(states,inputs,constants)>
%              (permitting '=', '==', or ':=' for assignment)
%    states (struct) - states of LTI system
%    inputs (struct) - input signals
%    parseMode - affects default outputs ('flow' (default), 'assignment')
%
% Outputs:
%    isLinear - equations successfully linearized? (logical)
%    A - state matrix of LTI system
%    B - input matrix of LTI system
%    c - constant offset of LTI system
%    formalEqs - equations formatted to paste into a matlab dynamics file
%    containsInput - boolean specifying if the assignment contains an input
%
% Example: 
%    [isLinear,A,B,c,formalEqs,containsInput] = ...
%       eq2linSys('x'' == v + xErr & v'' == -g + vErr',...
%       struct('name',{'x','v'}),struct('name',{'xErr','vErr'}),...
%       struct('name',{'g'},'value',{9.81}))
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

% set default value
parseMode = setDefaultValues({'flow'},varargin);

if ~strcmp(parseMode,'flow') && ~strcmp(parseMode,'assignment')
    throw(CORAerror('CORA:converterIssue',...
        'Wrong parse mode in conversion of flow equations.'));
end

% create symbolic variables for states
numStates = length(states);
containsInput = false;

stateNames = strings(numStates,1);
varnames = cell(numStates,1);
for i=1:numStates
    stateNames(i) = states(i).name;
    % add 'L'/'R' as placeholders for parentheses
    varnames{i} = strcat('xL',num2str(i),'R');
end
x = sym(varnames);

%create symbolic variables for inputs
numInputs = length(inputs);

inputNames = strings(numInputs,1);
varnames = cell(numInputs,1);
for i=1:numInputs
    inputNames(i) = inputs(i).name;
    % add 'L'/'R' as placeholders for parentheses
    varnames{i} = strcat('uL',num2str(i),'R');
end
u = sym(varnames);

% substitute states into equation
exprs = applySymMapping(exprs,stateNames,x);

% substitute inputs into equation
exprs = applySymMapping(exprs,inputNames,u);

% check whether any variables except states (x) and inputs (u) remain
isSubset = all(ismember(symvar(exprs),[x;u]));
if ~isSubset
    % print detailed error message: unknown variables
    diff = setdiff(symvar(exprs),[x;u]);
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
        if strcmp(states(i).name,exprNames{j})
            eqIdx = j;
            break;
        end
    end
    if eqIdx ~= 0
        f(i,1) = exprs(j);
    else
        if strcmp(parseMode,'assignment')
            % ommitting assignments allowed per SX specification
            f(i,1) = x(i);
        else % (default) mode 'flow'
            f(i,1) = sym(0);
            warning("No flow specified for state %s!\n" + ...
                "assuming %s' = 0...",states(i).name,states(i).name);
        end
    end
end

%computing linear dependencies
A_sym = jacobian(f,x);
B_sym = jacobian(f,u);

%computing constant components, by setting x = u = 0
f_const = subs(f,[x;u],zeros(numStates + numInputs,1));
%evaluate as column vector
c_sym = f_const(1:numStates);
    
try
    A = double(A_sym);
    B = double(B_sym);
    c = double(c_sym);
    isLinear = true;
    containsInput = any(any(B));
catch
    % Could not linearize system.
    A = []; B = []; c = [];
    isLinear = false;
    containsInput = any(ismember(symvar(exprs),u));
    
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
    % replace 'L' and 'R' with parentheses
    eq = strrep(eq,'L','(');
    eq = strrep(eq,'R',')');
    eqs(i) = sprintf("dx(%d,1) = %s;\n",i,eq);
end

formalEqs = join(eqs,newline);

% ------------------------------ END OF CODE ------------------------------
