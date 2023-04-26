function result = savepathonly(outputfile)
%SAVEPATHONLY Save the current MATLAB path to the file pathdef.m
%
%   SAVEPATHONLY saves the current MATLAB path to the current pathdef.m file.
%
%   SAVEPATHONLY(OUTPUTFILE) saves the current MATLAB path to the
%   file specified by OUTPUTFILE.
%
%   This is essentially a passthrough to SAVEPATH or PATH2RC, depending on
%   the Matlab version.
%
%   SAVEPATHONLY returns the result of SAVEPATH or PATH2RC:
%     0 if the file was saved successfully
%     1 if the file could not be saved
%
%   There are a few exceptional cases, which are not handled by SAVEPATH or
%   PATH2RC but should be.
%   1. The output of PATHDEF, MATLABPATH and PATH includes the directories
%   listed by USERPATH. SAVEPATH (or PATH2RC) does not check and remove
%   these directories. When the path is saved, SAVEPATH (or PATH2RC) adds
%   these extra directories (if any) to the pathdef.m file. The result when
%   Matlab is next started is that there are duplicate directories in the
%   path.
%   2. On Unix, if the directory "$HOME/matlab" is in the current MATLAB
%   path, then SAVEPATH (or PATH2RC) adds this directory to the pathdef.m
%   file, along with all the others in the Matlab path. This also results
%   in duplicate directories when Matlab is next started. This is because the
%   Matlab startup automatically adds the directory "$HOME/matlab" to the
%   environment variable MATLABPATH, if the directory exists. The function
%   USERPATH then uses the list provided by MATLABPATH, and as noted above,
%   the output of PATHDEF, etc. includes the directories listed by
%   USERPATH.
%
%   SAVEPATHONLY deals with these exceptional cases by removing the exceptional
%   directories from the path before calling SAVEPATH (or PATH2RC).
%
%   See also SAVEPATH (Matlab 7 or greater) or PATH2RC (Matlab 6.X).

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-05-31 $
% Initialize cell array removed{} as size 0 rather than size 1
% On Windows, do not remove [getenv('USERPROFILE') '\matlab'] directory
% $Revision 1.00 $ $Date 2005-02-03 $

%
% We will call call SAVEPATH (or PATH2RC) depending on Matlab version
%
matlab_ver = version;
if matlab_ver(1) >= '7'
    verge7 = true;
    savefn = @savepath;
else
    verge7 = false;
    savefn = @path2rc;
end
%
% Before version 7, only Unix has case sensitive path names
%
casesen = isunix || verge7;
p = path;
u = userpath;
if ~casesen
    p = lower(p);
    u = lower(u);
end

if ~nargin
    %
    % Using the current Matlab path, find pathdef.m
    %
    outputfile = which('pathdef');
end
%
% See if we need to also remove $HOME/matlab from the path
%
matlab_dir = '';
if isunix
    [status,dir_str] = system('echo $HOME');
    home_dir = deblank(dir_str);
    matlab_dir = fullfile(home_dir,'matlab');
end
if (~isempty(matlab_dir)) && isempty(strfind(u,[matlab_dir pathsep]))
    u = [matlab_dir pathsep u];
    if ~casesen
        u = lower(u);
    end
end
%
% Use the cell array REMOVED to keep track of directories rmoved
%
removed = cell(0);
%
% Remove directories from the path
%
k = 1;
while any(u)
    %
    % Get the next directory name from U and remove if necessary
    %
    [dirname u] = strtok(u, pathsep);
    if ~casesen
        dirname = lower(dirname);
    end
    if ~isempty(dirname) && ~isempty(strfind(p,[dirname pathsep]));
        % disp(['removing "' dirname '"']);
        rmpath(dirname);
        removed{k} = dirname;
        k = k+1;
    end
end
%
% Now that we have cleaned up the path we can safely call SAVEPATH (or PATH2RC)
%
result = feval(savefn,outputfile);
%
% We now add the paths back in
%
for k = 1:length(removed)
    % disp(['adding "' removed{k} '"']);
    addpath(removed{k});
end
