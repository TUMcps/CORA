function uninstall_eq_toolbox(arg)
%UNINSTALL_EQ_TOOLBOX Uninstall using Toolbox Installer.
%
%Syntax
% uninstall_eq_toolbox(arg);
%
%Description
% UNINSTALL_EQ_TOOLBOX uses Toolbox Installer to uninstall this toolbox.
% It must be run from the top level directory of this toolbox, which must
% contain the file info.ins.
%
% If UNINSTALL_EQ_TOOLBOX is called with no argument, it asks the user for
% confirmation before taking any action.
%
% UNINSTALL_EQ_TOOLBOX first checks if the file pathname.m exists on the current
% Matlab path and can be appended to.
% If so, it asks the user for confirmation and then simply calls the UNINSTALL
% function of Toolbox Installer to update the file pathname.m.
%
% UNINSTALL_EQ_TOOLBOX DIRNAME, that is UNINSTALL_EQ_TOOLBOX('DIRNAME'),
% calls UNINSTALL('DIRNAME'), using DIRNAME as the directory for pathdef.m.
%
% UNINSTALL_EQ_TOOLBOX -, that is UNINSTALL_EQ_TOOLBOX('-'), calls UNINSTALL('-'),
% which displays a dialog box so the user can select the directory to be used for
% pathdef.m.
%
%Notes
% UNINSTALL_EQ_TOOLBOX does not undo everything done by INSTALL_EQ_TOOLBOX.
% In particular, UNINSTALL_EQ_TOOLBOX does not delete pathdef.m or its
% directory, even if these were created by INSTALL_EQ_TOOLBOX.
%
% For Toolbox Installer 2.2 by B. Rasmus Anthin, see
% Matlab Central File Exchange
% http://www.mathworks.com/matlabcentral/fileexchange/
% Home page of B. Rasmus Anthin
% http://www.etek.chalmers.se/~e8rasmus/
%
%Examples
% > cd ~
% > uninstall_eq_toolbox
% ??? Error using ==> uninstall_eq_toolbox
% Please run uninstall_eq_toolbox from toolbox directory, eg. /c/leopardi/[...]
%
% > cd /c/leopardi/[...]
% > uninstall_eq_toolbox
% Info: Uninstaller will update the file /c/leopardi/matlab/pathdef.m
% Proceed (y/n)? y
% Info: Uninstalled.
%
% > uninstall_eq_toolbox
% Info: Uninstaller will update the file /c/leopardi/matlab/pathdef.m
% Proceed (y/n)? n
%
%See also
% INSTALL_EQ_TOOLBOX, PATHDEF, PATH, Toolbox Installer 2.2.

% Copyright 2004-2005 Paul Leopardi for the University of New South Wales.
% $Revision 1.10 $ $Date 2005-06-01 $
% Clearly distinguish between info and error messages
% Documentation files renamed
% $Revision 1.00 $ $Date 2005-02-13 $
%
% For licensing, see COPYING.
% For references, see AUTHORS.
% For revision history, see CHANGELOG.

command = mfilename;
[command_dir name ext ver] = fileparts(mfilename('fullpath'));
wd = pwd;
if ~strcmp(command_dir,wd)
    error(['Please run ' command ' from toolbox directory, eg. ' command_dir])
end

if nargin
    %
    % If user has used an argument, just call install with this argument.
    %
    install(arg);
    return;
end

estr = 'Error: ';
istr = 'Info: ';
was_not_uninstalled_msg = 'Uninstaller cannot continue.';
need_to_specify_msg     = 'You will need to specify the directory for pathdef.m';

pathdef_exists = (exist('pathdef') == 2);
if pathdef_exists
    pathdef_name = which('pathdef');
    pathdefid = fopen(pathdef_name,'a');
    if pathdefid > 0
        fclose(pathdefid);
        reply = input([istr 'Uninstaller will update the file ' strrep(pathdef_name,'\','\\') ...
            '\nProceed (y/n)? '],'s');
        if strcmpi(reply(1),'y')
            uninstall;
        end
        return;
    else
        disp([estr was_not_uninstalled_msg]);
        disp([estr 'Cannot write to the pathdef file ' ...
            pathdef_name ' which was found on the Matlab path.']);
        error(['Try "' command ' -". ' need_to_specify_msg]);
    end    
else
    disp([estr was_not_uninstalled_msg ' The file pathdef.m was not found.']);
    error(['Try "' command ' -". ' need_to_specify_msg]);
end
%
%end function
