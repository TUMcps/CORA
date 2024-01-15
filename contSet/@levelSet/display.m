function display(ls)
% display - Displays the properties of a levelSet object on the command
%    window
%
% Syntax:
%    display(ls)
%
% Inputs:
%    ls - levelSet object
%
% Outputs:
%    ---
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

% Authors:       Mark Wetzlinger
% Written:       09-June-2020
% Last update:   03-March-2022 (MP, levelSets with multiple equations)
%                15-May-2023 (MW, fix printed string)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if ~isempty(inputname(1))
    fprintf(newline);
    disp(inputname(1) + " =");
    fprintf(newline);
end

% string of variables in level set
varsPrintStr = strjoin(string(ls.vars),",");

% number of concatenated level sets
numSets = length(ls.eq);

if numSets == 1
    disp("  f(" + varsPrintStr + "): " ...
        + string(ls.eq) + " " + ls.compOp + " 0" + newline);
else
    % loop over all level sets
    for i=1:numSets
        % only use & sign for all but last equation
        andStr = "";
        if i ~= numSets
            andStr = " & ";
        end
        % display i-th equation
        disp("  f" + i + "(" + varsPrintStr + "): " ...
            + string(vpa(ls.eq(i),3)) + " " + ls.compOp(i) + " 0" ...
            + andStr);
    end
    fprintf(newline);
end

% ------------------------------ END OF CODE ------------------------------
