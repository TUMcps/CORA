function res = installCORA()
% installCORA - this script guides you through the installation process of
%    CORA
%
% Syntax:
%    installCORA
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Authors:       Tobias Ladner
% Written:       23-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init
res = false;

% add CORA repository to path
disp("Adding CORA to the MATLAB path..")
aux_add_corapath()

% display information
disp("--------------------------------------------------------------------------------")
% cora version
fprintf(['  Version:    <strong>',CORAVERSION,'</strong>\n']);
% matlab version
fprintf(['  Matlab:     ',version,'\n']);
% host system
hostsystem = 'Unknown';
if ismac
    hostsystem = 'Mac';
elseif isunix
    hostsystem = 'Unix';
elseif ispc
    hostsystem = 'Windows';
end
fprintf(['  System:     ',hostsystem,'\n']);
disp("--------------------------------------------------------------------------------")
disp(" ")

% install toolboxes
disp('Installing required toolboxes (see also <a href="cora.in.tum.de/manual">CORA manual</a>):')
disp('- [CORA] Installing required toolboxes..')

try 

% check yalmip
aux_install_yalmip()

% further toolboxes
aux_install_toolbox('Symbolic_Toolbox','Symbolic Math Toolbox','SM');
aux_install_toolbox('Optimization_Toolbox','Optimization Toolbox','OP');
aux_install_toolbox('Statistics_Toolbox','Statistics and Machine Learning Toolbox','ST');

% neural network verification
fprintf('- [CORA NNV] (Optional) Should the required toolboxes for the verification of neural networks be installed?\n')
fprintf('  (y,[n]) ')
answer = input('>> ', 's');
installnn = ~isempty(answer) && startsWith(answer,'y');
if installnn
    aux_install_toolbox('Neural_Network_Toolbox','Deep Learning Toolbox','NN');
    aux_install_supportpkg('Deep Learning Toolbox Converter for ONNX Model Format','ONNXCONVERTER');
end

catch ME
    disp(' ')
    if strcmp(ME.identifier,'CORA:install')
        % rethrow exception to remove call stack
        throw(CORAerror('CORA:install',ME.message));
    else
        rethrow(ME)
    end
end

% check if polytope is found
if ~startsWith(which('polytope'),CORAROOT)
    disp('<strong>CORA installation failed.</strong> CORA v2024 comes with a new polytope class and thus does not longer depend on the Multi Parametric Toolbox (MPT). Please remove MPT from the MATLAB path.')
    res = false;
    return;
end

% check installation
disp(' ')
if installnn
    res = testnn_requiredToolboxes;
else
    res = test_requiredToolboxes;
end
if res
    disp('<strong>CORA installed successfully!</strong> Please check the <a href="cora.in.tum.de/manual">CORA manual</a> and our <a href="cora.in.tum.de">website</a> to get started.')
else
    disp('<strong>CORA installation failed.</strong> Please rerun the script or check Sec. 1.3 in the <a href="cora.in.tum.de/manual">CORA manual</a> for help.')
end

end


% Auxiliary functions -----------------------------------------------------

function aux_add_corapath()
% get CORA root directory
coraroot = fileparts(which(mfilename));
addpath(genpath(coraroot));

end

function res = aux_install_toolbox(id_long,text,id_short)
    res = aux_test_toolbox_installation(id_long,text);
    while ~res
        aux_display_install_prompt(text,id_short);
    end
end

function res = aux_test_toolbox_installation(id,text)
    res = any(any(contains(struct2cell(ver),text)));
end

function aux_display_install_prompt(text, id_short)
    text = sprintf(['''<strong>%s</strong>'' is missing and requires manual installation. \n' ...
        ' 1. Please install it via the MATLAB <a href="matlab:matlab.internal.addons.launchers.showExplorer(''CORA'',''identifier'',''%s'')">Add-On Explorer</a>. \n' ...
        ' 2. <a href="matlab:installCORA">Rerun</a> the CORA installation script afterward.\n'], text, id_short);
    throw(CORAerror('CORA:install',text))
end

function aux_install_yalmip()
    if ~isYalmipInstalled()
        disp("  - Yalmip is installed via the MATLAB tbxmanager. Where should the tbxmanager be saved?")
        % partially taken from: https://www.mpt3.org/Main/Installation

        % find tbxmanager directory
        currentpath = pwd;
        tbxpath = uigetdir(currentpath);
        tbxpath = [tbxpath,filesep,'tbxmanager'];
        mkdir(tbxpath)
        cd(tbxpath);
        userpath(tbxpath);

        % from https://tbxmanager.com/
        websave('tbxmanager.m','http://www.tbxmanager.com/tbxmanager.m');
        rehash;

        errortext = [];

        try
            % install yalmip & sedumi
            tbxmanager install yalmip sedumi
            
        catch ME
            
            % for mac users, sedumi cannot be installed via the tbxmanager
            % trying to install only yalmip

            % install yalmip
            tbxmanager install yalmip

            % and tell user to install e.g. SDPT3
            errortext = ['Unfortunately, the sdpt solver ''sedumi'' could not be installed automatically together with yalmip.\n' ...
                'Please install another solver manually, e.g. SDPT3 https://blog.nus.edu.sg/mattohkc/softwares/sdpt3/'];

        end

        % modify startup
        startuppath = which('startup.m');
        if isempty(startuppath)
            startuppath = [tbxpath,filesep,'startup.m'];
        end
        fid = fopen(startuppath,'a');
        if isequal(fid,-1)
            error(['Could not modify the initialization file "startup.m".',...
                   'Edit this file in the folder "%s" manually and insert there the line:  tbxmanager restorepath.'],startuppath);
        end

        % write startup file
        fprintf(fid,'%% startup\n');
        fprintf(fid,'fprintf(''Startup file: %%s\\n'',which(mfilename))\n');
        fprintf(fid,'\n');
        fprintf(fid,'%% init tbxmanager\n');
        fprintf(fid,'tbxmanager restorepath\n');
        fprintf(fid,'\n');
        fprintf(fid,'%% init cora\n');
        fprintf(fid,'addpath(genpath(''%s''));\n', CORAROOT);
        fprintf(fid,'fprintf(''Toolbox "%%s" added to the Matlab path.\\n'',CORAVERSION)\n');
        fprintf(fid,'\n');
        fprintf(fid,'disp(''Done.'')\n');
        fclose(fid);

        % change directory back to current directory
        cd(currentpath);
        addpath(tbxpath);

        % throw error message if sedumi is missing
        if ~isempty(errortext)
            throw(CORAerror('CORA:install',errortext));
        end
    end
end

function res = aux_install_supportpkg(text,id_short)
    res = aux_test_supportpkg_installation(text);
    while ~res
        aux_display_install_prompt(text,id_short);
    end
end

function res = aux_test_supportpkg_installation(text)
    addons = matlabshared.supportpkg.getInstalled;
    res = ~isempty(addons) && any(strcmp(text, [addons.Name]));
end

% ------------------------------ END OF CODE ------------------------------
