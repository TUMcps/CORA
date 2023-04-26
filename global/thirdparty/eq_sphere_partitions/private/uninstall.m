function uninstall(pddir)
%UNINSTALL  Uninstall toolbox directories.
%   UNINSTALL is an easy to use uninstallation program for toolboxes.
%   It uninstalls toolbox directories when run from the base directory
%   of a toolbox. A file "info.ins" created by MAKEINSTALL must exist
%   in this directory.
%   Example of usage:
%
%       cd(fullfile(matlabroot,'toolbox','digitalsim'))
%       disp(pwd)
%      C:\MATLAB\toolbox\digitalsim
%       uninstall
%      Uninstalled.
%      
%
%   If you are sharing matlab over a network, the location of the
%   'pathdef.m' file might have to be specified:
%   UNINSTALL PDDIR
%   In the example above, we could instead write:
%
%       uninstall C:\MATLAB\toolbox\local\
%
%   Alternatively, you can write:
%   UNINSTALL -
%   which will bring up an file selection menu where you choose the
%   pathdef.m file you want UNINSTALL to update.
%
%   See also INSTALL, MAKEINSTALL, CHECKINSTALL, SAVEPATHONLY.

% Revision 2005-05-28 Copyright (c) Paul Leopardi for UNSW:
% Clearly distinguish between info and error messages
% Revision 2005-02-01 Copyright (c) Paul Leopardi for UNSW:
% Use rmpath only when directory name is in path.
% Do not warn on IS_INSTALLED; do not save info.ins
% Use lower case if not case sensitive; case sesitive if Unix or version >= 7
% Use savepathonly rather than path2rc
% Display message if we did not need to remove any directories from the path
% Copyright (c) 2003-07-15, B. Rasmus Anthin.
% Revision 2003-07-16, 2003-07-17, 2003-07-22.

istr = 'Info: ';
if ~exist(fullfile(pwd,'info.ins'),'file')
   error('This toolbox cannot be uninstalled.')
end
load info.ins -mat
%% if ~IS_INSTALLED
%%    warning('This toolbox is already uninstalled.')
%% end

matlab_ver = version;
verge7 = (matlab_ver(1) >= '7');
%
% Before version 7, only Unix has case sensitive path names
%
casesen = isunix || verge7;

p = path;
wd = pwd;
if ~casesen
    wd = lower(wd);
    p = lower(p);
end
removed_any = false;
for i = 1:length(INS_DIRS)
    dirname = fullfile(wd,INS_DIRS{i});
    if ~casesen
        dirname = lower(dirname);
    end
    if ~isempty(strfind(p,[dirname pathsep]));
        rmpath(dirname);
        removed_any = true;
    end
end
if ~isempty(strfind(p,[wd pathsep]))
    rmpath(wd);
    removed_any = true;
end

if ~removed_any
    disp([istr 'Uninstaller did not need to remove any directories from the Matlab path.']);
end

if ~nargin
   flag=savepathonly;
elseif strcmp(pddir,'-')
   [fname,pname]=uigetfile('pathdef.m','Choose ''pathdef.m'' to update.');
   flag=savepathonly(fullfile(pname,fname));
else
   flag=savepathonly(fullfile(pddir,'pathdef.m'));
end
switch(flag)
case 1, error('Path could not be saved to pathdef.m. Try "uninstall -".')
case 2, error('Original pathdef.m was not found. Try "uninstall -".')
case 3, error('Original pathdef.m was found but couldn''t be read. Try "uinstall -".')
end
disp([istr 'Uninstalled.'])
%% IS_INSTALLED=0;
%% save info.ins INS_DIRS IS_INSTALLED
