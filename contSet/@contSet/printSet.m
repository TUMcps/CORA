function res = printSet(varargin)
% printSet - prints a set such that if one executes this command
%    in the workspace, this set would be created
%
% Syntax:
%    printSet(S)
%    printSet(S,accuracy,doCompact,clearLine)
%    printSet(fid,S,varargin)
%    printSet(filename,S,varargin)
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
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contDynamics/printSystem, printMatrix, printCell, printStruct

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   20-May-2025 (TL, added option for fid)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,S,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% get print info
[abbrev,propertyOrder] = getPrintSetInfo(S);

if doCompact
    % print in one line
    fprintf(fid,'%s(',class(S));
    for p = 1:numel(propertyOrder)
        pname = propertyOrder{p};
        aux_printProperty(fid,S.(pname),accuracy);
        if p < numel(propertyOrder)
            fprintf(fid,', ');
        end
    end
    fprintf(fid,')');

else
    % print each property as variable
    for p = 1:numel(propertyOrder)
        pname = propertyOrder{p};
        fprintf(fid,'%s = ',pname);
        aux_printProperty(fid,S.(pname),accuracy);
        fprintf(fid,';\n');
    end
    % init set
    fprintf(fid,'%s = %s(%s);',abbrev,class(S),strjoin(propertyOrder,','));
end

if clearLine
    fprintf(fid,'\n');
end

% finalize
res = closePrint(fid,closefid);

end


% Auxiliary functions -----------------------------------------------------

function aux_printProperty(fid,property,accuracy)
    if isnumeric(property)
        printMatrix(fid,property,accuracy,true,false);
    elseif iscell(property)
        printCell(fid,property,accuracy,true,false);
    elseif isstruct(property)
        printStruct(fid,property,accuracy,true,false);
    elseif isa(property,'contSet')
        printSet(fid,property,accuracy,true,false)
    else
        throw(CORAerror("CORA:noops",property))
    end
end


% ------------------------------ END OF CODE ------------------------------
