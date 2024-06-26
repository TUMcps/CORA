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
%    clearLine - (optional) whether to finish with '\n'
%    doCompact - (optional) whether the matrix is printed compactly
%
% Outputs:
%    -
%
% Example: 
%    M = [2 3; -2 1];
%    printMatrix(M)

% Authors:       Matthias Althoff
% Written:       01-November-2017
% Last update:   27-June-2018
%                04-January-2021
%                17-June-2022 (MW, parsing of accuracy)
%                06-June-2024 (TL, added doCompact)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin > 4
    throw(CORAerror('CORA:tooManyInputArgs',4));
end

% determine accuracy
if nargin < 2 || isempty(varargin{1})
    % default accuracy
    accuracy = '%4.3f%s';
else
    % numerical value
    if isnumeric(varargin{1})
        accuracy = varargin{1};
    % category
    elseif ischar(varargin{1})
        if strcmp(varargin{1},'high')
            accuracy = '%16.16f%s';
        else
            accuracy = varargin{1};
        end
    end
end

% determine clearLine
if nargin < 3
    clearLine = true;
else
    clearLine = varargin{2};
end

if nargin < 4
    doCompact = false;
else
    doCompact = varargin{3};
end

% scalar case
if isscalar(M)
    fprintf(accuracy, M);
    if clearLine
        fprintf('\n')
    end
    return
end

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
        %write in workspace
        fprintf(accuracy, M(iRow,iCol));

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
