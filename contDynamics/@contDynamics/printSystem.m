function res = printSystem(varargin)
% printSystem - prints a contDynamic object such that if one executes this command
%    in the workspace, this contDynamic object would be created
%
% Syntax:
%    printSystem(sys)
%    printSystem(sys,accuracy,doCompact,clearLine)
%    printSystem(fid,sys,varargin)
%    printSystem(filename,sys,varargin)
%
% Inputs:
%    sys - contDynamics
%    accuracy - floating-point precision
%    doCompact - whether to compactly print the set
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
% See also: printSet, printMatrix, printCell, printStruct

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   20-May-2025 (TL, added option for fid)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,sys,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% get print info
[propertyOrder] = getPrintSystemInfo(sys);

if doCompact
    % print in one line
    fprintf(fid,'%s(',class(sys));
    for p = 1:numel(propertyOrder)
        pname = propertyOrder{p};
        aux_printProperty(fid,sys.(pname),accuracy);
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
        aux_printProperty(fid,sys.(pname),accuracy);
        fprintf(fid,';\n');
    end
    % init set
    fprintf(fid,'sys = %s(%s);',class(sys),strjoin(propertyOrder,','));
end

if clearLine
    fprintf(fid,'\n');
end

% finalize
res = closePrint(fid,closefid);

end


% Auxiliary functions -----------------------------------------------------

function aux_printProperty(fid,property,accuracy)
    % Matlab objects ---
    if ischar(property) || isstring(property)
        fprintf(fid,"'%s'",property);
    elseif isnumeric(property)
        printMatrix(fid,property,accuracy,true,false);
    elseif iscell(property)
        printCell(fid,property,accuracy,true,false);
    elseif isstruct(property)
        printStruct(fid,property,accuracy,true,false);
    elseif isa(property,'function_handle')
        funcstr = func2str(property);
        if ~startsWith(funcstr,'@')
            funcstr = sprintf('@%s',funcstr);
        end
        fprintf(fid,funcstr);
        % CORA objects ---
    elseif isa(property,'contSet') || isa(property,'matrixSet')
        printSet(fid,property,accuracy,true,false);
    elseif isa(property,'contDynamics')
        printSystem(fid,property,accuracy,true,false);
    else
        throw(CORAerror("CORA:noops",property))
    end
end


% ------------------------------ END OF CODE ------------------------------
