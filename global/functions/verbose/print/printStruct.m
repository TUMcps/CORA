function res = printStruct(varargin)
% printStruct - prints a struct such that if one executes this command
%    in the workspace, this struct would be created
%
% Syntax:
%    printStruct(S)
%    printStruct(S,accuracy,doCompact,clearLine)
%    printStruct(fid,S,varargin)
%    printStruct(filename,S,varargin)
%
% Inputs:
%    S - struct array
%    accuracy - floating-point precision
%    doCompact - whether to compact matrices
%    clearLine - whether to finish with '\n'
%    filename - char, filename to print given obj to
%    fid - char, fid to print given obj to
%
% Outputs:
%    res - logical
%
% Example: 
%    S = struct('a',[1 2 3],'b','text');
%    printStruct(S)
%
% See also: printMatrix, printCell, printSet

% Authors:       Tobias Ladner
% Written:       31-May-2023
% Last update:   14-May-2024 (TL, added doCompact)
%                10-October-2024 (TL, clean up input parsing)
%                20-May-2025 (TL, added option for fid)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,S,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% find num rows and cols
numRows = size(S,1);
numCols = size(S,2);

% write first element
if numRows > 1 
    fprintf(fid,'[ ');
    if ~doCompact
        fprintf(fid,'...\n ');
    end
end

% iterate over struct array
for i = 1:numRows
    for j = 1:numCols
        % read current struct
        S_ij = S(i,j);
        fprintf(fid,'struct(');

        % iterate over fieldnames
        names = fieldnames(S_ij);
        numNames = length(names);

        if numNames > 1 && ~doCompact
            fprintf(fid,' ...\n');
        end

        for k = 1:numNames
            name= names{k};
            if numNames > 1 && ~doCompact
                % intend
                fprintf(fid,'    ');
            end
            fprintf(fid,"'%s', ", name);

            % print value s.t. it can be re-created
            value = S_ij.(name);
            if isnumeric(value) || islogical(value)
                printMatrix(fid,value,accuracy,doCompact,false);
            elseif isstruct(value)
                printStruct(fid,value,accuracy,doCompact,false);
            elseif iscell(value)
                % enclosing with {...} due to struct constructor
                fprintf(fid,'{');
                printCell(fid,value,accuracy,doCompact,false);
                fprintf(fid,'}');
            elseif ischar(value) || isstring(value)
                fprintf(fid,"'%s'",value);
            elseif isa(value,'contSet')
                % always print compactly
                printSet(fid,value,accuracy,true,false);
            elseif isa(value,'contDynamics')
                % always print compactly
                printSystem(fid,value,accuracy,true,false);
            else
                % print as string and hope for the best
                fprintf(fid,"'%s'",value);
            end
            if k<numNames
                fprintf(fid,',');
            end

            if numNames > 1 && ~doCompact
                fprintf(fid,' ...\n');
            end
        end
        fprintf(fid,')');
        if j < numCols
            fprintf(fid,',');
        end
    end
    if length(S) > 1
        if doCompact
            fprintf(fid,'; ');
        else
            fprintf(fid,'; ...\n');
        end
    end
end

% finish struct
if length(S) > 1
    fprintf(fid,']');
end

if clearLine
    fprintf(fid,'\n');
end

% finalize
res = closePrint(fid,closefid);

end


% ------------------------------ END OF CODE ------------------------------
