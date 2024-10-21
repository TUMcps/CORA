function printStruct(S,varargin)
% printStruct - prints a struct such that if one executes this command
%    in the workspace, this struct would be created
%
% Syntax:
%    printStruct(S)
%    printStruct(S,'high')
%
% Inputs:
%    S - struct array
%    accuracy - (optional) floating-point precision
%    doCompact - (optional) whether to compact matrices
%    clearLine - (optional) whether to finish with '\n'
%
% Outputs:
%    -
%
% Example: 
%    S = struct('a',[1 2 3],'b','text');
%    printStruct(S)
%
% See also: printMatrix, printCell, printSet

% Authors:       Tobias Ladner
% Written:       31-May-2023
% Last update:   14-May-2024 (TL, added coCompact)
%                10-October-2024 (TL, clean up input parsing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(0,4);
[accuracy,doCompact,clearLine] = setDefaultValues({'%4.3f%s',false,true},varargin);
if ischar(accuracy) && strcmp(accuracy,'high')
    accuracy = '%16.16f%s';
end
inputArgsCheck({ ...
    {S,'att','struct'}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

numRows = size(S,1);
numCols = size(S,2);

% write first element
if length(S) > 1
    fprintf('[ ...\n ')
end

% iterate over struct array
for i = 1:numRows
    for j = 1:numCols
        % read current struct
        S_ij = S(i,j);
        fprintf('struct(')

        % iterate over fieldnames
        names = fieldnames(S_ij);
        numNames = length(names);

        if numNames > 1 && ~doCompact
            fprintf(' ...\n')
        end

        for k = 1:numNames
            name= names{k};
            if numNames > 1 && ~doCompact
                % intend
                fprintf('    ')
            end
            fprintf("'%s', ", name)

            % print value s.t. it can be re-created
            value = S_ij.(name);
            if isnumeric(value)
                printMatrix(value,accuracy,doCompact,false)
            elseif isstruct(value)
                printStruct(value,accuracy,doCompact,false)
            elseif iscell(value)
                fprintf('{')
                printCell(value,accuracy,doCompact,false)
                fprintf('}')
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
            if k<numNames
                fprintf(',')
            end

            if numNames > 1 && ~doCompact
                fprintf(' ...\n')
            end
        end
        fprintf(')')
        if j < numCols
            fprintf(',')
        end
    end
    if length(S) > 1 && ~doCompact
        fprintf('; ...\n')
    end
end

% finish struct
if length(S) > 1
    fprintf(']')
end

if clearLine
    fprintf('\n')
end


% ------------------------------ END OF CODE ------------------------------
