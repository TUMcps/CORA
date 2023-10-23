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
%    clearLine - (optional) whether to finish with '\n'
%
% Outputs:
%    -
%
% Example: 
%    S = struct('a',[1 2 3],'b','text');
%    printStruct(S)
%

% Authors:       Tobias Ladner
% Written:       31-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% determine accuracy
if nargin < 2
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

        if numNames > 1
            fprintf(' ...\n')
        end

        for k = 1:numNames
            name= names{k};
            if numNames > 1
                % intend
                fprintf('    ')
            end
            fprintf("'%s', ", name)

            % print value s.t. it can be re-created
            value = S_ij.(name);
            if isnumeric(value)
                printMatrix(value,accuracy,false)
            elseif isa(value,'interval')
                printInterval(value)
            elseif isa(value,"struct")
                printStruct(value,accuracy,false)
            elseif ischar(value) || isstring(value)
                fprintf("'%s'",value)
            else
                % print as string and hope for the best
                fprintf("%s",value)
            end
            if k<numNames
                fprintf(',')
            end

            if numNames > 1
                fprintf(' ...\n')
            end
        end
        fprintf(')')
        if j < numCols
            fprintf(',')
        end
    end
    if length(S) > 1
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
