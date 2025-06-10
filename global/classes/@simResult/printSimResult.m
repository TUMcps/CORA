function res = printSimResult(varargin)
% printSimResult - prints a simResult object such that if one executes this command
%    in the workspace, this simResult object would be created
%
% Syntax:
%    printSimResult(simRes)
%    printSimResult(simRes,accuracy,doCompact,clearLine)
%    printSimResult(fid,simRes,varargin)
%    printSimResult(filename,simRes,varargin)
%
% Inputs:
%    simRes - simResult object
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
% See also: specification

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   20-May-2025 (TL, added option for fid)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,simRes,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% call constructor
fprintf(fid,'simResult(');
if ~doCompact
    fprintf(fid,' ...\n');
end
% x
printCell(fid,simRes.x,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% t
printCell(fid,simRes.t,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% loc
loc = simRes.loc;
if iscell(loc)
    loc = cell2mat(loc);
end
printMatrix(fid,loc,'%i',true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% y
printCell(fid,simRes.y,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% a
printCell(fid,simRes.a,accuracy,true,false);
if ~doCompact
    fprintf(fid,' ...\n');
end
fprintf(fid,')');

if clearLine
    fprintf(fid,'\n');
end

% finalize
res = closePrint(fid,closefid);

end

% ------------------------------ END OF CODE ------------------------------
