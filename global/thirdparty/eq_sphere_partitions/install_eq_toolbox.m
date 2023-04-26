function install_eq_toolbox(arg)
%INSTALL_EQ_TOOLBOX Install using Toolbox Installer, with sensible defaults
%
%Syntax
% install_eq_toolbox(arg);
%
%Description
% INSTALL_EQ_TOOLBOX uses Toolbox Installer to install this toolbox. It must
% be run from the top level directory of this toolbox, which must contain the
% file info.ins.
%
% If INSTALL_EQ_TOOLBOX is called with no argument, it asks the user for
% confirmation before taking any action.
%
% INSTALL_EQ_TOOLBOX first checks if the file pathname.m exists on the current
% Matlab path and can be appended to.
% If so, it asks the user for confirmation and then simply calls the INSTALL
% function of Toolbox Installer to update the file pathname.m.
% If not, it tries to choose a sensible directory in which to create a new
% pathdef.m file, and suggests this to the user before proceeding.
%
% INSTALL_EQ_TOOLBOX DIRNAME, that is INSTALL_EQ_TOOLBOX('DIRNAME'),
% calls INSTALL('DIRNAME'), using DIRNAME as the directory for pathdef.m.
%
% INSTALL_EQ_TOOLBOX -, that is INSTALL_EQ_TOOLBOX('-'), calls INSTALL('-'),
% which displays a dialog box so the user can select the directory to be used
% for pathdef.m.
%
%Notes
% If it is necessary to create a new pathdef.m file, the directory which
% INSTALL_EQ_TOOLBOX suggests to the user for this file defaults to
% the directory which contains startup.m, if this exists.
%
% Otherwise INSTALL_EQ_TOOLBOX suggests the following directory to the user:
% For Unix systems including Linux and Mac: $HOME/matlab
% For Windows systems: [matlabroot '\toolbox\local'].
% If this directory does not exist, INSTALL_EQ_TOOLBOX tries to create it.
%
% For Toolbox Installer 2.2 by B. Rasmus Anthin, see
% Matlab Central File Exchange
% http://www.mathworks.com/matlabcentral/fileexchange/
% Home page of B. Rasmus Anthin
% http://www.etek.chalmers.se/~e8rasmus/
%
%Examples
% > install_eq_toolbox
% Info: Installer will create the file /c/leopardi/matlab/pathdef.m
% Proceed (y/n)? y
% Info: Installer has created the directory /c/leopardi/matlab
% Info: Installed.
%
% > install_eq_toolbox
% Info: Installer will update the file /c/leopardi/matlab/pathdef.m
% Proceed (y/n)? y
% Info: Installed.
%
% > install_eq_toolbox
% Info: Installer will update the file /c/leopardi/matlab/pathdef.m
% Info: Proceed (y/n)? n
%
% > cd ~
% > install_eq_toolbox
% ??? Error using ==> install_eq_toolbox
% Please run install_eq_toolbox from toolbox directory, eg. /c/leopardi/[...]
%
%See also
% UNINSTALL_EQ_TOOLBOX, PATHDEF, PATH, Toolbox Installer 2.2

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Partial overhaul
% o Restructure tests for existing pathdef.m file and guesses for pathdef_dir
% o Simplify tests after mkdir
% o On Windows, do not use [getenv('USERPROFILE') '\matlab'] directory
% o Clearly distinguish between info and error messages
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-12 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

if (nargin == 1) && strcmp(arg,'-')
    install(arg);
    return;
end

command = mfilename;
[command_dir name ext ver] = fileparts(mfilename('fullpath'));
wd = pwd;
if ~strcmp(command_dir,wd)
    error(['Please run ' command ' from toolbox directory, eg. ' command_dir]);
end

if nargin
    %
    % If user has used an argument, just call install with this argument.
    %
    install(arg);
    return;
end

estr = 'Error: ';
wstr = 'Warning: ';
istr = 'Info: ';
was_not_installed_msg = 'Installer cannot continue.';
try_to_create_msg =     'Try to create the directory outside of Matlab';

current_pathdef_name = which('pathdef');
pathdef_exists = (exist('pathdef') == 2);
if pathdef_exists
    pathdef_name = current_pathdef_name;
    %
    % We have been told that pathdef.m exists.
    % Try to open it for updating
    %
    pathdefid = fopen(pathdef_name,'a');
    if pathdefid > 0
        %
        % File open for pathdef.m succeeded.
        % Ask whether to proceed with update.
        %
        fclose(pathdefid);
        reply = input([istr 'Installer will update the file ' strrep(pathdef_name,'\','\\') ...
            '\nProceed (y/n)? '],'s');
        if strcmpi(reply(1),'y')
            install;
        end
        return;
    elseif (exist('startup') == 2)
        %
        % File open for pathdef.m failed.
        % We have been told that startup.m exists.
        % Guess that we should use pathdef.m in the same directory.
        %
        startup_name = which('startup');
        [pathdef_dir name ext ver] = fileparts(startup_name);
        pathdef_name = fullfile(pathdef_dir,'pathdef.m');
    else
        %
        % File open for pathdef.m failed.
        % We have been told that no startup.m exists.
        %
        if isunix
            %
            % We are using Unix.
            % Guess that we should use $HOME/matlab/pathdef.m
            %
            [status,dir_str] = system('echo $HOME');
            home_dir = deblank(dir_str);
            pathdef_name = fullfile(home_dir,'matlab','pathdef.m');
        else
            %
            % We are using Windows.
            % Guess that we should use $matlabroot\toolbox\local\pathdef.m
            %
            pathdef_name = fullfile(matlabroot,'toolbox','local','pathdef.m');
        end
    end
    [pathdef_dir name ext ver] = fileparts(pathdef_name);
end

if strcmp(pathdef_name, current_pathdef_name)
    %
    % If we have not changed pathdef_name,
    % then this is the same file that we have been told we cannot write to,
    % and we have not been able to guess another possiblity for pathdef_name.
    %
    disp([wstr 'Cannot write to the pathdef file ' ...
        pathdef_name ' which was found on the Matlab path.']);
    reply = input([istr 'Do you want to select an existing pathdef.m file to update (y/n)? '], 's');
    if strcmpi(reply(1),'y')
        %
        % The user has chosen to update an existing pathdef.m file.
        % Use install -
        %
        install('-');
        return;
    end
    reply = input([istr 'Do you want to select a directory and create a new pathdef.m file (y/n)? '], 's');
    if ~strcmpi(reply(1),'y')
        %
        % The user has chosen not to create a new pathdef.m file.
        % There are no other possibilities, so return.
        %
        return;
    end
    %
    % The user has chosen to create a new pathdef.m file.
    % Ask where it should go. The directory must exist.
    %
    pname=uigetdir('','Please choose an existing directory for the new ''pathdef.m'' file.');
    if isequal(pname,0)
        %
        % The user has cancelled the directory selection.
        % There are no other possibilities, so return.
        %
        return;
    end
    pathdef_name = fullfile(pname,'pathdef.m');
    [pathdef_dir name ext ver] = fileparts(pathdef_name);
end
%
% Determine if Matlab thinks that pathdef_name already exists as a file.
%
if exist(pathdef_name) == 2
    willstr = 'update';
else
    willstr = 'create';
end
%
% Ask before attempting to perform the appropriate action.
%
reply = input([istr 'Installer will ' willstr ' the file ' strrep(pathdef_name,'\','\\') ...
    '\nProceed (y/n)? '],'s');
if ~strcmpi(reply(1),'y')
    return;
end
%
% Determine if Matlab thinks that pathdef_dir already exists as a directory.
%
pathdef_dir_exists = (exist(pathdef_dir) == 7);
%
% Even if Matlab says that pathdef_dir exists, it might not exist as a directory.
% Try to create it anyway. This should not hurt if it really does exist.
%
[parent dirname ext ver] = fileparts(pathdef_dir);
[status message] = mkdir(parent,dirname);
if status < 0
    created_pathdef_dir = false;
    if ~pathdef_dir_exists
        %
        % Matlab said pathdef_dir does not exist and that we cannot not create it.
        %
        disp([estr was_not_installed_msg]);
        disp([estr 'Cannot create ' pathdef_dir ': ' message]);
        error([try_to_create_msg '\n and then try "' command pathdef_dir '"']);
    end
else
    created_pathdef_dir = true;
    if ~pathdef_dir_exists
        %
        % Matlab has said that pathdef_dir did not exist and that we have now created it.
        %
        disp([istr command ' has created the directory ' pathdef_dir]);
    end
end
%
% At this point, pathdef_dir should exist as a directory.
% Check again.
%
pathdef_dir_exists = (exist(pathdef_dir) == 7);
if ~pathdef_dir_exists
    %
    % Matlab said pathdef_dir does not exist.
    %
    disp([estr was_not_installed_msg]);
    disp([estr 'Cannot create ' pathdef_dir ': ' message]);
    error([try_to_create_msg '\n and then try "' command pathdef_dir '"']);
end
%
% Finally, we can ask install to use pathdef_dir.
%
install(pathdef_dir);
%
% If we created pathdef_dir, we should add it to the path,
% so that the current Matlab session can see the new pathdef.m
%
if created_pathdef_dir
    addpath(pathdef_dir);
end
rehash pathreset;
%
%end function
