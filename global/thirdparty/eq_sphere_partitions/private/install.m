function install(pddir)
%INSTALL  Install toolbox directories.
%   INSTALL is an easy to use installation program for toolboxes.
%   It installs toolbox directories when run from the base directory
%   of a toolbox. A file "info.ins" created by MAKEINSTALL must exist
%   in this directory.
%   Example of usage:
%
%       cd(fullfile(matlabroot,'toolbox','digitalsim'))
%       disp(pwd)
%      C:\MATLAB\toolbox\digitalsim
%       install
%      Installed.
%
%
%   If you are sharing matlab over a network, the location of the
%   'pathdef.m' file might have to be specified:
%   INSTALL PDDIR
%   In the example above, we could instead write:
%
%       install C:\MATLAB\toolbox\local\
%
%   Alternatively, you can write:
%   INSTALL -
%   which will bring up an file selection menu where you choose the
%   pathdef.m file you want INSTALL to update.
%
%   See also UNINSTALL, MAKEINSTALL, CHECKINSTALL, SAVEPATHONLY.

% Revision 2005-05-28 Copyright (c) Paul Leopardi for UNSW:
% Clearly distinguish between info and error messages
% Check return value of uigetfile and return if 0.
% Revision 2005-02-01 Copyright (c) Paul Leopardi for UNSW:
% Do not warn on IS_INSTALLED; do not save info.ins
% Use savepathonly rather than path2rc.
% Copyright (c) 2003-07-15, B. Rasmus Anthin.
% Revision 2003-07-16, 2003-07-17, 2003-07-22.

istr = 'Info: ';
if ~exist(fullfile(pwd,'info.ins'),'file')
   error('This toolbox cannot be installed.')
end
load info.ins -mat
%% if IS_INSTALLED
%%    warning('This toolbox is already installed.')
%% end
addpath(pwd)
for i=1:length(INS_DIRS)
   addpath(fullfile(pwd,INS_DIRS{i}))
end
if ~nargin
   flag=savepathonly;
elseif strcmp(pddir,'-')
   [fname,pname]=uigetfile('pathdef.m','Choose ''pathdef.m'' to update.');
   if isequal(fname,0) | isequal(pname,0)
      return;
   end
   flag=savepathonly(fullfile(pname,fname));
else
   flag=savepathonly(fullfile(pddir,'pathdef.m'));
end
switch(flag)
case 1, error('Path could not be saved to pathdef.m. Try "install -".')
case 2, error('Original pathdef.m was not found. Try "install -".')
case 3, error('Original pathdef.m was found but couldn''t be read. Try "install -".')
end
disp([istr 'Installed.'])
%% IS_INSTALLED=1;
%% save info.ins INS_DIRS IS_INSTALLED
