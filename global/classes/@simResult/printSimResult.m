function printSimResult(simRes,varargin)
% printSimResult - prints a simResult object such that if one executes this command
%    in the workspace, this simResult object would be created
%
% Syntax:
%    printSimResult(simResult)
%    printSimResult(simResult,'high')
%
% Inputs:
%    simRes - simResult object
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
    {simRes,'att','simResult'}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

fprintf('simResult(')
if ~doCompact
    fprintf(' ...\n')
end
printCell(simRes.x,accuracy,true,false)
fprintf(', ')
if ~doCompact
    fprintf('...\n')
end
printCell(simRes.t,accuracy,true,false)
fprintf(', ')
if ~doCompact
    fprintf('...\n')
end
loc = simRes.loc;
if iscell(loc)
    loc = cell2mat(loc);
end
printMatrix(loc,'%i',true,false)
fprintf(', ')
if ~doCompact
    fprintf('...\n')
end
printCell(simRes.y,accuracy,true,false)
fprintf(', ')
if ~doCompact
    fprintf('...\n')
end
printCell(simRes.a,accuracy,true,false)
if ~doCompact
    fprintf(' ...\n')
end
fprintf(')')

if clearLine
    fprintf('\n')
end

end

% ------------------------------ END OF CODE ------------------------------
