function [ A, b, Ae, be ] = eq2polytope(IneqExprs, EqExprs, states)
% Converts string inequations to a polytope specification (if possible).
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
 %      A - matrix of size: #inequations X length(states)
 %      b - vector of size: length(states)
 %      Ae - matrix of size: #equations X length(states)
 %      be - vector of size: length(states)
 %      Represented polytope: A*x <= b
 %      Represented polyhedron: A*x <= b & Ae*x == b
 %
 % EXAMPLE: eq2polytope('x <= eps & v < 0', struct('name',{'x','v'}),struct('name',{'eps'},'value',{'0.75'}))
 
 %create symbolic variables for states
 numStates = length(states);
 
 stateNames = strings(numStates,1);
 varnames = cell(numStates,1);
 for i=1:numStates
     stateNames(i) = states(i).name;
     varnames{i} = strcat('xL',num2str(i),'R');
 end
 x = sym(varnames);
 
 % substitute state variables into equations
 IneqExprs = applySymMapping(IneqExprs,stateNames,x);
 EqExprs = applySymMapping(EqExprs,stateNames,x);
 
  % check wether any expressions include non-state variables
 % (compute a logical index for this)
 ineqValidIdx = false(size(IneqExprs));
 for i=1:length(ineqValidIdx)
     ineqValidIdx(i) = all(ismember(symvar(IneqExprs(i)), x));
 end
 eqValidIdx = false(size(EqExprs));
 for i=1:length(eqValidIdx)
     eqValidIdx(i) = all(ismember(symvar(EqExprs(i)), x));
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
         if ~eqValidIdx(i)
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
 
 try
    A = double(A_sym);
    b = double(b_sym) * -1; %reminder: IneqExprs = A*x - b
    
    Ae = double(Ae_sym);
    be = double(be_sym) * -1; %reminder: EqExprs = Ae*x - be
 catch
     error('Could not express as Polyhedron (not linear):\n%s',equation);
 end
end