function res = printCell(varargin)
% printCell - prints a cell array such that if one executes this command
%    in the workspace, this cell array would be created
%
% Syntax:
%    printCell(C)
%    printCell(C,accuracy,doCompact,clearLine)
%    printCell(fid,C,varargin)
%    printCell(filename,C,varargin)
%
% Inputs:
%    C - cell array
%    accuracy - floating-point precision
%    doCompact - whether the matrix is printed compactly
%    clearLine - whether to finish with '\n'
%    filename - char, filename to print given obj to
%    fid - char, fid to print given obj to
%
% Outputs:
%    res - logical
%
% Example: 
%    C = {2 3; -2 1};
%    printCell(C)
%
% See also: printMatrix, printStruct, printSet

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   20-May-2025 (TL, added option for fid)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,C,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% empty case
if isempty(C)
    if all(size(C) == 0)
        fprintf(fid,'{}');
    else
        fprintf(fid,'cell(%i,%i)',size(C,1),size(C,2));
    end
    if clearLine
        fprintf(fid,'\n');
    end
    res = closePrint(fid,closefid);
    return
end

% general case
numRows = size(C,1);
numCols = size(C,2);

% write first element
fprintf(fid,'{ ');
if numRows > 1 && ~doCompact
    fprintf(fid,'...\n ');
end

%write each row
for iRow=1:numRows
    for iCol=1:numCols
        %write in command window
        value = C{iRow,iCol};
        if isnumeric(value) || islogical(value)
            printMatrix(fid,value,accuracy,doCompact,false);
        elseif isa(value,"struct")
            printStruct(fid,value,accuracy,doCompact,false);
        elseif ischar(value) || isstring(value)
            fprintf(fid,"'%s'",value);
        elseif isa(value,'contSet')
            % always print compactly (as no longer one-liner)
            printSet(fid,value,accuracy,true,false);
        elseif isa(value,'contDynamics')
            % always print compactly (as no longer one-liner)
            printSystem(fid,value,accuracy,true,false);
        elseif isa(value,'simResult')
            printSimResult(fid,value,accuracy,doCompact,false);
        else
            % print as string and hope for the best
            fprintf(fid,"%s",value);
        end

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
fprintf(fid,'}');
if clearLine
    fprintf(fid,'\n');
end
res = closePrint(fid,closefid);

end

% ------------------------------ END OF CODE ------------------------------
