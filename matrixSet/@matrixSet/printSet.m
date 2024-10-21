function printSet(S,varargin)
% printSet - prints a set such that if one executes this command
%    in the workspace, this set would be created
%
% Syntax:
%    printSet(S)
%    printSet(S,'high')
%
% Inputs:
%    S - contSet
%    accuracy - (optional) floating-point precision
%    doCompact - (optional) whether to compactly print the set
%    clearLine - (optional) whether to finish with '\n'
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/printSystem, printMatrix, printCell, printStruct

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

narginchk(0,4);
[accuracy,doCompact,clearLine] = setDefaultValues({'%4.3f%s',false,true},varargin);
if ischar(accuracy) && strcmp(accuracy,'high')
    accuracy = '%16.16f%s';
end
inputArgsCheck({ ...
    {S,'att',{'contSet','matrixSet'}}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

% get print info
[abbrev,propertyOrder] = getPrintSetInfo(S);

if doCompact
    % print in one line
    fprintf('%s(',class(S))
    for p = 1:numel(propertyOrder)
        pname = propertyOrder{p};
        aux_printProperty(S.(pname),accuracy);
        if p < numel(propertyOrder)
            fprintf(', ')
        end
    end
    fprintf(')')

else
    % print each property as variable
    for p = 1:numel(propertyOrder)
        pname = propertyOrder{p};
        fprintf('%s = ',pname);
        aux_printProperty(S.(pname),accuracy);
        fprintf(';\n');
    end
    % init set
    fprintf('%s = %s(%s);',abbrev,class(S),strjoin(propertyOrder,','));
end

if clearLine
    fprintf('\n');
end

end


% Auxiliary functions -----------------------------------------------------

function aux_printProperty(property,accuracy)
    if isnumeric(property)
        printMatrix(property,accuracy,true,false);
    elseif iscell(property)
        printCell(property,accuracy,true,false);
    elseif isstruct(property)
        printStruct(property,accuracy,true,false);
    elseif isa(property,'contSet')
        printSet(property,accuracy,true,false)
    else
        throw(CORAerror("CORA:noops",property))
    end
end


% ------------------------------ END OF CODE ------------------------------
