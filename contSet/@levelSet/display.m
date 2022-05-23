function display(ls)
% display - Displays the levelSet object
%
% Syntax:  
%    display(ls)
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    (to console)
%
% Example: 
%    syms x y
%    eq = x^2 + y^2 - 4;
%    ls = levelSet(eq,[x;y],'==');
%    display(ls);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      09-June-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

fprintf(newline);
disp(inputname(1) + " =");
fprintf(newline);

varsStr = string(ls.vars);
nrVars = length(varsStr);
varsPrintStr = "";
for i=1:nrVars-1
    varsPrintStr = varsPrintStr + varsStr(i) + ",";
end
varsPrintStr = varsPrintStr + varsStr(nrVars);

if length(ls.compOp) == 1
    disp("  eq: f(" + varsPrintStr + ") " + ls.compOp + " " + string(ls.eq) + newline);
else
    for i = 1:length(ls.compOp)
        andStr = "";
        if i ~= length(ls.compOp)
            andStr = " & ";
        end
        disp("  eq_"+i+": f(" + varsPrintStr + ") " + ls.compOp(i) + " " + string(vpa(ls.eq(i),3)) + andStr + newline);
    end
end

%------------- END OF CODE --------------