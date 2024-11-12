function printCell(C,varargin)
% printCell - prints a cell array such that if one executes this command
%    in the workspace, this cell array would be created
%
% Syntax:
%    printCell(C)
%
% Inputs:
%    C - cell array
%    accuracy - (optional) floating-point precision
%    doCompact - (optional) whether the matrix is printed compactly
%    clearLine - (optional) whether to finish with '\n'
%
% Outputs:
%    -
%
% Example: 
%    C = {2 3; -2 1};
%    printCell(C)
%
% See also: printMatrix, printStruct, printSet

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(0,4);
[accuracy,doCompact,clearLine] = setDefaultValues({'%4.3f',false,true},varargin);
if ischar(accuracy) && strcmp(accuracy,'high')
    accuracy = '%16.16f';
end
inputArgsCheck({ ...
    {C,'att','cell'}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

% empty case
if isempty(C)
    if all(size(C) == 0)
        fprintf('{}');
    else
        fprintf('cell(%i,%i)',size(C,1),size(C,2));
    end
    if clearLine
        fprintf('\n')
    end
    return
end

% general case
numRows = size(C,1);
numCols = size(C,2);

% write first element
fprintf('{ ');
if numRows > 1 && ~doCompact
    fprintf('...\n ')
end

%write each row
for iRow=1:numRows
    for iCol=1:numCols
        %write in command window
        value = C{iRow,iCol};
        if isnumeric(value)
            printMatrix(value,accuracy,doCompact,false)
        elseif isa(value,"struct")
            printStruct(value,accuracy,doCompact,false)
        elseif ischar(value) || isstring(value)
            fprintf("'%s'",value)
        elseif isa(value,'contSet')
            % always print compactly
            printSet(value,accuracy,true,false)
        elseif isa(value,'contDynamics')
            % always print compactly
            printSystem(value,accuracy,true,false)
        else
            % print as string and hope for the best
            fprintf("%s",value)
        end

        if iCol < numCols && ~doCompact
            fprintf(',')
        end
        fprintf(' ')
    end

    if numRows > 1 || iRow<numRows
        if doCompact
            if iRow < numRows
                fprintf('; ')
            end
        else
            %write new line
            fprintf('; ...\n ')
        end
    end
end
fprintf('}')
if clearLine
    fprintf('\n')
end

% ------------------------------ END OF CODE ------------------------------
