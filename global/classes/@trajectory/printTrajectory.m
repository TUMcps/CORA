function res = printTrajectory(varargin)
% printTrajectory - prints a trajectory object such that if one executes this command
%    in the workspace, this trajectory object would be created
%
% Syntax:
%    printTrajectory(traj)
%    printTrajectory(traj,accuracy,doCompact,clearLine)
%    printTrajectory(fid,traj,varargin)
%    printTrajectory(filename,traj,varargin)
%
% Inputs:
%    traj - trajectory object
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

% Authors:       Tobias Ladner, Laura Luetzow
% Written:       08-August-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
[fid,closefid,traj,accuracy,doCompact,clearLine] = initPrint(varargin{:});

% call constructor
fprintf(fid,'trajectory(');
if ~doCompact
    fprintf(fid,' ...\n');
end
% u
printMatrix(fid,traj.u,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% x
printMatrix(fid,traj.x,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% y
printMatrix(fid,traj.y,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% t
printMatrix(fid,traj.t,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% dt
printMatrix(fid,traj.dt,accuracy,true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% loc
loc = traj.loc;
if iscell(loc)
    loc = cell2mat(loc);
end
printMatrix(fid,loc,'%i',true,false);
fprintf(fid,', ');
if ~doCompact
    fprintf(fid,'...\n');
end
% a
printMatrix(fid,traj.a,accuracy,true,false);
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
