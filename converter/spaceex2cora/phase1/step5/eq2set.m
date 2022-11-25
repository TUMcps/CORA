function set = eq2set(IneqExprs, EqExprs, states, outputs)
% Converts string of inequalities to a CORA set representation
% Requires names of state variables.
% Requires names of constants and values to substitute.
% WARNING: could produce incorrect results, if any variables are named 
% "xL<number>R"

 % INPUT
 %      equation: string - rough format: <eq> ("&" <eq>)*
 %                  <eq> = <expr(states,constants)> <op> <expr(states,constants)>
 %                  <expr>: linear expressions only
 %                  <op> = "<"|">"|"<="|">="|==
 % OUTPUT
 %      set - mptPolytope or levelSet object
 %
 % EXAMPLE: eq2polytope('x <= eps & v < 0', struct('name',{'x','v'}),struct('name',{'eps'},'value',{'0.75'}))
 
 %create symbolic variables for states + outputs
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
 IneqExprs = applySymMapping(IneqExprs,listOfVarNames,x);
 EqExprs = applySymMapping(EqExprs,listOfVarNames,x);
 
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
 ineqValidIdx = false(size(IneqExprs));
 for i=1:length(ineqValidIdx)
     ineqValidIdx(i) = all(ismember(symvar(IneqExprs(i)), x_y));
 end
 eqValidIdx = false(size(EqExprs));
 for i=1:length(eqValidIdx)
     eqValidIdx(i) = all(ismember(symvar(EqExprs(i)), x_y));
 end
 
 % if any expressions fail the test, print warnings then remove them
 if(~all(ineqValidIdx) || ~all(eqValidIdx))
     
     %print warnings for inequalities
     for i=1:length(ineqValidIdx)
         if ~ineqValidIdx(i)
             % detect, which variables caused the error
             diff = setdiff(symvar(IneqExprs(i)), x);
             %print detailed warning message
             varstr = "(" + string(diff(1));
             if(length(diff) > 1)
                 varstr = varstr + sprintf(", %s",diff(2:end));
             end
             varstr = varstr + ")";
             warning("A CONDITION CONTAINS NON-STATE VARIABLES %s AND IS IGNORED!"+...
                    newline + "Condition: ""%s <= 0"" (after arithmetic transformation)",...
                    varstr,string(IneqExprs(i)));
         end
     end
     
     %print warnings for equalities
     for i=1:length(eqValidIdx)
         if ~eqValidIdx(i) && ~all(ismember(symvar(applySymMapping(EqExprs(i),listOfVarNames,x_y)), x_y))
             % detect, which variables caused the error
             diff = setdiff(symvar(EqExprs(i)), x);
             %print detailed warning message
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
     
     %remove invalid expressions from arrays
     IneqExprs = IneqExprs(ineqValidIdx);
     EqExprs = EqExprs(eqValidIdx);
 end
  
 %computing linear dependencies
 A_sym = jacobian(IneqExprs,x);
 Ae_sym = jacobian(EqExprs,x);
 
 %computing constant components, by setting x = 0
 b_sym = subs(IneqExprs , x, zeros(length(x),1));
 be_sym = subs(EqExprs , x, zeros(length(x),1));
 
 % equations are linear -> construct mptPolytope
 try
    Opts.A = double(A_sym);
    Opts.b = double(b_sym) * -1; %reminder: IneqExprs = A*x - b
    
    Opts.Ae = double(Ae_sym);
    Opts.be = double(be_sym) * -1; %reminder: EqExprs = Ae*x - be
    
    set = mptPolytope(Opts);
    
 % equations are nonlinear -> construct levelSet
 catch
     if ~isempty(EqExprs)
        if ~isempty(EqExprs)
            eq = [EqExprs;IneqExprs];
        else
            eq = IneqExprs; 
        end
     else
        eq = IneqExprs; 
     end
     
     compOp1 = repmat({'=='},[length(EqExprs),1]);
     compOp2 = repmat({'<='},[length(IneqExprs),1]);
     compOp = [compOp1;compOp2];
     
     if length(compOp) == 1
        set = levelSet(eq,x,compOp{1}); 
     else
        set = levelSet(eq,x,compOp); 
     end    
 end
end