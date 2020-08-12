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

disp("  eq: f(" + varsPrintStr + ") " + ls.compOp + " " + string(ls.eq) + newline);


%------------- END OF CODE --------------