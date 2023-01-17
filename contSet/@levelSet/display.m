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

% Author:       Mark Wetzlinger
% Written:      09-June-2020
% Last update:  03-March-2022 (MP, add functionality for levelSets with
%                                 multiple equations)
% Last revision:---

%------------- BEGIN CODE --------------

if isemptyobject(ls)
    
    dispEmptyObj(ls,inputname(1));

else

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
        disp("  eq: f(" + varsPrintStr + ") " + ls.compOp + " " ...
            + string(ls.eq) + newline);
    else
        % loop over all level sets
        for i=1:numSets
            % only use & sign for all but last equation
            andStr = "";
            if i ~= numSets
                andStr = " & ";
            end
            % display i-th equation
            disp("  eq. #" + i + ": f(" + varsPrintStr + ") " + ls.compOp(i) ...
                + " " + string(vpa(ls.eq(i),3)) + andStr);
        end
    end

end

%------------- END OF CODE --------------