function printSpec(spec,varargin)
% printSpec - prints a specification object such that if one executes this command
%    in the workspace, this specification object would be created
%
% Syntax:
%    printSpec(spec)
%    printSpec(spec,'high')
%
% Inputs:
%    spec - specification
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
% See also: specification

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
    {spec,'att','specification'}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

for i=1:numel(spec)
    fprintf('specification(')
    if ~doCompact
        fprintf(' ...\n')
    end
    printSet(spec(i).set,accuracy,true,false)
    fprintf(', ')
    if ~doCompact
        fprintf('...\n')
    end
    fprintf("'%s'",spec(i).type)
    fprintf(', ')
    if ~doCompact
        fprintf('...\n')
    end
    printSet(spec(i).time,accuracy,true,false)
    fprintf(', ')
    if ~doCompact
        fprintf('...\n')
    end
    printMatrix(spec(i).location,'%i',true,false)
    if ~doCompact
        fprintf(' ...\n')
    end
    fprintf(')')

    if clearLine
        fprintf('\n')
    end
end

end

% ------------------------------ END OF CODE ------------------------------
