function res = printMatrix(varargin)
% printMatrix - prints an matrix such that if one executes this command
%    in the workspace, this matrix would be created
%
% Syntax:
%    printMatrix(M)
%    printMatrix(M,accuracy,doCompact,clearLine)
%    printMatrix(fid,M,varargin)
%    printMatrix(filename,M,varargin)
%
% Inputs:
%    M - matrix
%    accuracy - floating-point precision, or 'high'
%    doCompact - whether the matrix is printed compactly
%    clearLine - whether to finish with '\n'
%    filename - char, filename to print given obj to
%    fid - char, fid to print given obj to
%
% Outputs:
%    res - logical
%
% Example: 
%    M = [2 3; -2 1];
%    printMatrix(M)
%
% See also: printCell, printStruct, printSet

% Authors:       Matthias Althoff, Tobias Ladner
% Written:       01-November-2017
% Last update:   27-June-2018
%                04-January-2021
%                17-June-2022 (MW, parsing of accuracy)
%                06-June-2024 (TL, added doCompact)
%                10-October-2024 (TL, clean up input parsing)
%                20-May-2025 (TL, added option for fid,optimizations)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,M,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% convert logical values to numeric
M = M * 1;

% check if all integer
if isequal(M,round(M))
    accuracy = '%i';
end

% empty case
if isempty(M)
    if all(size(M) == 0)
        fprintf(fid,'[]');
    else
        fprintf(fid,'zeros(%i,%i)',size(M,1),size(M,2));
    end
    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% scalar case
if isscalar(M)
    if round(M,0) == M % check integer
        fprintf(fid,'%i', M);
    else
        fprintf(fid,accuracy, M);
    end
    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% matrix of one value case
if isscalar(full(unique(M)))
    dims = size(M);
    dimstr = sprintf(join(repmat("%i",1,numel(dims)),','),dims);
    value = full(M(1));
    if value == 0
        mat_string = sprintf('zeros(%s)',dimstr);
    elseif value == 1
        mat_string = sprintf('ones(%s)',dimstr);
    else
        % display using value*ones(n,m)
        format_string = strcat(accuracy,"*ones(%s)");
        mat_string = sprintf(format_string,full(M(1,1)),dimstr);
    end
    if issparse(M)
        fprintf(fid,"sparse(%s)\n",mat_string);
    else
        fprintf(fid,mat_string);
    end
    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% sparse case
if issparse(M)
    [i,j,v] = find(M);
    fprintf(fid,'sparse(');
    if ~doCompact
        fprintf(fid,' ...\n');
    end
    printMatrix(i','%i',true,false);
    fprintf(fid,', ');
    if ~doCompact
        fprintf(fid,'...\n');
    end
    printMatrix(j','%i',true,false);
    fprintf(fid,', ');
    if ~doCompact
        fprintf(fid,'...\n');
    end
    printMatrix(v',accuracy,true,false);
    if ~doCompact
        fprintf(fid,' ...\n');
    end
    fprintf(fid,')');

    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% nd case
if numel(size(M)) > 2
    % reshape to 1d vector and print with reshape to nd array
    fprintf(fid,'reshape(');
    printMatrix(fid,reshape(M,1,[]),accuracy,true,false);
    % compact size in reshape call
    fprintf(fid,', [%s])',strrep(num2str(size(M)),'   ',','));
    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% 2d case ---

numRows = size(M,1);
numCols = size(M,2);

% check if (partial) identity
if isequal(M,eye(numRows,numCols))
    if numCols == numRows
        mat_string = sprintf('eye(%i)',numRows);
    else
        mat_string = sprintf('eye(%i,%i)',numRows,numCols);
    end
    fprintf(fid,mat_string);
    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% write first element
fprintf(fid,'[ ');
if numRows > 1 && ~doCompact
    fprintf(fid,'...\n ');
end

%write each row
for iRow=1:numRows
    for iCol=1:numCols
        %write in command window
        fprintf(fid,accuracy, full(M(iRow,iCol)));

        if iCol < numCols && ~doCompact
            fprintf(fid,',');
        end
        fprintf(fid,' ');
    end

    if numRows > 1 || iRow<numRows
        if doCompact
            if iRow < numRows
                fprintf(fid,'; ');
            end
        else
            %write new line
            fprintf(fid,'; ...\n ');
        end
    end
end
fprintf(fid,']');
if clearLine
    fprintf(fid,'\n');
end
res = closePrint(fid,closefid);

end

% ------------------------------ END OF CODE ------------------------------
