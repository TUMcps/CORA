function printMatrix(M,varargin)
% printMatrix - prints an matrix such that if one executes this command
%    in the workspace, this matrix would be created
%
% Syntax:
%    printMatrix(M)
%
% Inputs:
%    M - matrix
%    accuracy - (optional) floating-point precision
%    doCompact - (optional) whether the matrix is printed compactly
%    clearLine - (optional) whether to finish with '\n'
%
% Outputs:
%    -
%
% Example: 
%    M = [2 3; -2 1];
%    printMatrix(M)
%
% See also: printCell, printStruct, printSet

% Authors:       Matthias Althoff
% Written:       01-November-2017
% Last update:   27-June-2018
%                04-January-2021
%                17-June-2022 (MW, parsing of accuracy)
%                06-June-2024 (TL, added doCompact)
%                10-October-2024 (TL, clean up input parsing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(0,4);
[accuracy,doCompact,clearLine] = setDefaultValues({'%4.3f',false,true},varargin);
if ischar(accuracy) && strcmp(accuracy,'high')
    accuracy = '%16.16f';
end
inputArgsCheck({ ...
    {M,'att','numeric'}, ...
    {accuracy,'att', {'char','string'}}, ...
    {doCompact,'att','logical'}, ...
    {clearLine,'att','logical'}
})

% empty case
if isempty(M)
    if all(size(M) == 0)
        fprintf('[]');
    else
        fprintf('zeros(%i,%i)',size(M,1),size(M,2));
    end
    if clearLine
        fprintf('\n')
    end
    return
end

% scalar case
if isscalar(M)
    if round(M,0) == M % check integer
        fprintf('%i', M);
    else
        fprintf(accuracy, M);
    end
    if clearLine
        fprintf('\n')
    end
    return
end

% matrix of one value case
if isscalar(full(unique(M)))
    [n_rows,n_clmns] = size(M);
    if M(1,1) == 0
        mat_string = sprintf('zeros(%i,%i)',n_rows,n_clmns);
    else
        % use user-specified accuracy
        format_string = strcat(accuracy,"*ones(%i,%i)");
        mat_string = sprintf(format_string,full(M(1,1)),n_rows,n_clmns);
    end
    if issparse(M)
        fprintf("sparse(%s)\n",mat_string);
    else
        fprintf(mat_string);
    end
    if clearLine
        fprintf('\n');
    end
    return
end

% sparse case
if issparse(M)
    [i,j,v] = find(M);
    fprintf('sparse(')
    if ~doCompact
        fprintf(' ...\n')
    end
    printMatrix(i','%i',true,false);
    fprintf(', ')
    if ~doCompact
        fprintf('...\n')
    end
    printMatrix(j','%i',true,false);
    fprintf(', ')
    if ~doCompact
        fprintf('...\n')
    end
    printMatrix(v',accuracy,true,false);
    if ~doCompact
        fprintf(' ...\n')
    end
    fprintf(')')

    if clearLine
        fprintf('\n')
    end
    return
end

% nd case
if numel(size(M)) > 2
    % reshape to 1d vector and print with reshape to nd array
    fprintf('reshape(')
    printMatrix(reshape(M,1,[]),accuracy,true,false);
    % compact size in reshape call
    fprintf(', [%s])',strrep(num2str(size(M)),'  ',','))
    if clearLine
        fprintf('\n')
    end
    return
end

% 2d case
numRows = size(M,1);
numCols = size(M,2);

% write first element
fprintf('[ ');
if numRows > 1 && ~doCompact
    fprintf('...\n ')
end

%write each row
for iRow=1:numRows
    for iCol=1:numCols
        %write in command window
        fprintf(accuracy, full(M(iRow,iCol)));

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
fprintf(']')
if clearLine
    fprintf('\n')
end

% ------------------------------ END OF CODE ------------------------------
