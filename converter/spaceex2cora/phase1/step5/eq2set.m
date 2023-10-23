function set = eq2set(ineqExprs,eqExprs,states,outputs)
% eq2set - Converts string of inequalities to a CORA set representation
%    Requires names of state variables.
%    Requires names of constants and values to substitute.
%    WARNING: could produce incorrect results, if any variables are named 
%             "xL<number>R"
%
% Syntax:
%    set = eq2set(ineqExprs,eqExprs,states,outputs)
%
% Inputs:
%    ineqExprs - symbolic expressions containing inequalities 
%                (of the form 'term <= 0')
%    eqExprs - symbolic expressions containing equalities
%              (of the form 'term == 0')
%    states - list of states of the component which the equation belongs to
%    outputs - list of outputs of the component which the equation belongs to
%
% Outputs:
%    set - polytope or levelSet object
%
% Example: 
%    set = eq2set('x <= eps & v < 0',struct('name',{'x','v'}),...
%       struct('name',{'eps'},'value',{'0.75'}));
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

% create symbolic variables for states + outputs
numStates = length(states);
if ~isempty(outputs)
    numOutputs = length(outputs);
else
    numOutputs = 0;
end

listOfVarNames = strings(numStates,1);
varNames = cell(numStates,1);
for i=1:numStates
    listOfVarNames(i) = states(i).name;
    varNames{i} = strcat('xL',num2str(i),'R');
end
x = sym(varNames);

% substitute state variables into equations
ineqExprs = applySymMapping(ineqExprs,listOfVarNames,x);
eqExprs = applySymMapping(eqExprs,listOfVarNames,x);

% include outputs to varNames, otherwise error if
%  output specified in invariant (has to be done so in spaceex...)
j = i+1;
for i=1:numOutputs
    listOfVarNames(j,1) = outputs(i).name;
    varNames{j,1} = strcat('yL',num2str(i),'R');
    j = j+1;
end
x_y = sym(varNames);


% check whether any expressions include non-state variables
% (compute a logical index for this)
ineqValidIdx = false(size(ineqExprs));
for i=1:length(ineqValidIdx)
    ineqValidIdx(i) = all(ismember(symvar(ineqExprs(i)), x_y));
end
eqValidIdx = false(size(eqExprs));
for i=1:length(eqValidIdx)
    eqValidIdx(i) = all(ismember(symvar(eqExprs(i)), x_y));
end

% if any expressions fail the test, print warnings then remove them
if(~all(ineqValidIdx) || ~all(eqValidIdx))
 
    %print warnings for inequalities
    for i=1:length(ineqValidIdx)
        if ~ineqValidIdx(i)
            % detect, which variables caused the error
            diff = setdiff(symvar(ineqExprs(i)), x);
            %print detailed warning message
            varstr = "(" + string(diff(1));
            if(length(diff) > 1)
                varstr = varstr + sprintf(", %s",diff(2:end));
            end
            varstr = varstr + ")";
            warning("A CONDITION CONTAINS NON-STATE VARIABLES %s AND IS IGNORED!"+...
            newline + "Condition: ""%s <= 0"" (after arithmetic transformation)",...
            varstr,string(ineqExprs(i)));
        end
    end
     
    %print warnings for equalities
    for i=1:length(eqValidIdx)
        if ~eqValidIdx(i) && ~all(ismember(symvar(applySymMapping(eqExprs(i),listOfVarNames,x_y)), x_y))
            % detect, which variables caused the error
            diff = setdiff(symvar(eqExprs(i)), x);
            %print detailed warning message
            varstr = "(" + string(diff(1));
            if(length(diff) > 1)
                varstr = varstr + sprintf(", %s",diff(2:end));
            end
            varstr = varstr + ")";
            warning("A CONDITION CONTAINS NON-STATE VARIABLES %s AND IS IGNORED!"+...
            newline + "Condition: ""%s == 0"" (after arithmetic transformation)",...
            varstr,string(eqExprs(i)));
        end
    end
    
    %remove invalid expressions from arrays
    ineqExprs = ineqExprs(ineqValidIdx);
    eqExprs = eqExprs(eqValidIdx);
end

%computing linear dependencies
A_sym = jacobian(ineqExprs,x);
Ae_sym = jacobian(eqExprs,x);

%computing constant components, by setting x = 0
b_sym = subs(ineqExprs,x,zeros(length(x),1));
be_sym = subs(eqExprs,x,zeros(length(x),1));

% equations are linear -> construct polytope
try
    % reminder: IneqExprs = A*x - b
    P_A = double(A_sym);
    P_b = double(b_sym) * -1;
    
    % reminder: EqExprs = Ae*x - be
    P_Ae = double(Ae_sym);
    P_be = double(be_sym) * -1;
    
    % instantiate polytope
    set = polytope(P_A,P_b,P_Ae,P_be);

% equations are nonlinear -> construct levelSet
catch
    if ~isempty(eqExprs)
        if ~isempty(eqExprs)
            eq = [eqExprs;ineqExprs];
        else
            eq = ineqExprs; 
        end
    else
        eq = ineqExprs; 
    end
    
    compOp1 = repmat({'=='},[length(eqExprs),1]);
    compOp2 = repmat({'<='},[length(ineqExprs),1]);
    compOp = [compOp1;compOp2];
    
    if length(compOp) == 1
        set = levelSet(eq,x,compOp{1}); 
    else
        set = levelSet(eq,x,compOp); 
    end    
end

% ------------------------------ END OF CODE ------------------------------
