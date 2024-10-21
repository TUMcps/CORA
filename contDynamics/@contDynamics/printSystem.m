function printSystem(sys,varargin)
% printSystem - prints a contDynamic object such that if one executes this command
%    in the workspace, this contDynamic object would be created
%
% Syntax:
%    printSystem(sys)
%    printSystem(sys,'high')
%
% Inputs:
%    sys - contDynamic
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
% See also: printSet, printMatrix, printCell, printStruct

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
    {sys,'att',{'contDynamics'}}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

% get print info
[propertyOrder] = getPrintSystemInfo(sys);

if doCompact
    % print in one line
    fprintf('%s(',class(sys))
    for p = 1:numel(propertyOrder)
        pname = propertyOrder{p};
        aux_printProperty(sys.(pname),accuracy);
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
        aux_printProperty(sys.(pname),accuracy);
        fprintf(';\n');
    end
    % init set
    fprintf('sys = %s(%s);',class(sys),strjoin(propertyOrder,','));
end

if clearLine
    fprintf('\n');
end

end


% Auxiliary functions -----------------------------------------------------

function aux_printProperty(property,accuracy)
    if ischar(property) || isstring(property)
        fprintf("'%s'",property);
    elseif isnumeric(property)
        printMatrix(property,accuracy,true,false);
    elseif iscell(property)
        printCell(property,accuracy,true,false);
    elseif isstruct(property)
        printStruct(property,accuracy,true,false);
    elseif isa(property,'contSet') || isa(property,'matrixSet')
        printSet(property,accuracy,true,false)
    elseif isa(property,'contDynamics')
        printSystem(property,accuracy,true,false)
    elseif isa(property,'function_handle')
        fprintf(func2str(property))
    else
        throw(CORAerror("CORA:noops",property))
    end
end


% ------------------------------ END OF CODE ------------------------------
