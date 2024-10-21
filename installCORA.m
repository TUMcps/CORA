function res = installCORA(varargin)
% installCORA - this script guides you through the installation process of
%    CORA
%
% Syntax:
%    installCORA
%
% Inputs:
%    interactive - logical, whether installation should be done interactively
%    installNNV - logical, whether toolboxed for NNV should be installed
%    defaultpath - char, path where toolboxes should be installed
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
% Last update:   10-January-2024 (added non-interactive mode)
%                29-January-2024 (changes for automatic installation in docker)
%                04-August-2024 (limit number of cpu threads)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% add CORA repository to path
disp("Adding CORA to the Matlab path..")
aux_add_corapath()

% init
[interactive,installNNV,defaultpath] = setDefaultValues({true,false,[]},varargin);

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
disp('  [1/4] YALMIP')
aux_install_yalmip(interactive,defaultpath)

% further toolboxes
toolboxes = {
    % id_long, title, id_short
    'Symbolic_Toolbox','Symbolic Math Toolbox','SM';
    'Optimization_Toolbox','Optimization Toolbox','OP';
    'Statistics_Toolbox','Statistics and Machine Learning Toolbox','ST';
};
for i=1:size(toolboxes,1)
    fprintf('  [%i/%i] %s\n', i+1,size(toolboxes,1)+1, toolboxes{i,2})
    aux_install_toolbox(toolboxes{i,:});
end

% neural network verification
if interactive
    fprintf('  (Optional) Should the required toolboxes for the verification of neural networks be installed?\n')
    fprintf('  (y,[n]) ')
    answer = input('>> ', 's');
    installNNV = ~isempty(answer) && startsWith(lower(answer),'y');
end
if installNNV
    disp('- [CORA NNV] Installing required toolboxes..')
    fprintf('  [1/3] %s\n', 'Deep Learning Toolbox')
    aux_install_toolbox('Neural_Network_Toolbox','Deep Learning Toolbox','NN');
    fprintf('  [2/3] %s\n', 'Deep Learning Toolbox Converter for ONNX Model Format')
    aux_install_supportpkg('Deep Learning Toolbox Converter for ONNX Model Format','ONNXCONVERTER');
    fprintf('  [3/3] %s\n', 'Parallel Computing Toolbox')
    aux_install_toolbox('Parallel_Computing_Toolbox','Parallel Computing Toolbox','DM');
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
    disp('<strong>CORA installation failed.</strong> CORA v2024 comes with a new polytope class and thus does not longer depend on the Multi Parametric Toolbox (MPT). Please remove MPT from the Matlab path.')
    res = false;
    return;
end

% check installation
disp(' ')
if installNNV
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
        ' 1. Please install it via the Matlab <a href="matlab:matlab.internal.addons.launchers.showExplorer(''CORA'',''identifier'',''%s'')">Add-On Explorer</a>. \n' ...
        ' 2. <a href="matlab:installCORA">Rerun</a> the CORA installation script afterward.\n'], text, id_short);
    throw(CORAerror('CORA:install',text))
end

function aux_install_yalmip(interactive,defaultpath)
    if ~isYalmipInstalled()
        disp("  - Yalmip is installed via the Matlab tbxmanager.")
        % partially taken from: https://www.mpt3.org/Main/Installation

        % find tbxmanager directory
        currentpath = pwd;
        if interactive && isempty(defaultpath)
            disp("    Where should the tbxmanager be saved?")
            disp("    (opening folder selection dialog box..)")
            pause(1)
            tbxpath = uigetdir(currentpath,'[CORA] Where should the tbxmanager be saved?');
        else
            tbxpath = defaultpath;
        end
        tbxpath = [tbxpath,filesep,'tbxmanager'];
        fprintf("    Saving tbxmanager to: %s\n",tbxpath);
        mkdir(tbxpath)
        cd(tbxpath);

        % from https://tbxmanager.com/
        if interactive
            % take original tbxmanger
            websave('tbxmanager.m','http://www.tbxmanager.com/tbxmanager.m');
        else
            % take tbxmanger from verivital to avoid manual license agreement
            websave('tbxmanager.m','https://raw.githubusercontent.com/verivital/tbxmanager/master/tbxmanager.m');
        end
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

        % also install mp
        % mppath = [fileparts(tbxpath),filesep,'MultiPrecisionToolbox'];
        % try
        %     fprintf('\nInstalling "MultiPrecision" toolbox to: %s\n\n',mppath);
        %     % download mp toolbox
        %     mkdir(mppath)
        %     mpzip = [mppath filesep 'mp.zip'];
        %     websave(mpzip,'https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/e54ba25f-4a80-11e4-9553-005056977bd0/4b93a1ca-e163-425c-a70e-345faa48d6bc/packages/zip');
        %     % unzip
        %     unzip(mpzip);
        %     % delete zip
        %     delete(mpzip);
        % catch ME
        %     fprintf('    MultiPrecision toolbox could not be installed automatically: %s', ME.message)
        % end

        % modify startup file
        startuppaths = {
            which('startup.m'); % append to existing
            [defaultpath filesep 'startup.m']; % default diretory
            [tbxpath filesep 'startup.m']; % append to tbxmanager
            [currentpath filesep 'startup.m']; % current diretory
        };
        for i=1:numel(startuppaths)
            startuppath = startuppaths{i};
            fid = fopen(startuppath,'a');
            if ~isequal(fid,-1)
                % able to write in this file
                break
            end
        end
        if isequal(fid,-1)
            % no place to write startup.m was found
            error('Could not modify the initialization file "startup.m".');
        end
        fprintf('Writing startup file to: %s\n', startuppath)

        % get startup lines
        startuptext = fileread(startuppath);
        lines = splitlines(startuptext);

        % adapt existing startup files
        isFunction = false;
        if ~isempty(lines) 
            % check if startup file starts with "function"
            for i=1:numel(lines)
                if startsWith(lines{i},'function startup')
                    % remove "end" at the end of function
                    isFunction = true;
                    break
                end
            end

            % remove "end"
            if isFunction
                for i=0:(length(lines)-1)
                    if contains(lines{end-i},'end')
                        lines{end-i} = '';
                        break
                    end
                end
            end
        end

        % write CORA startup file
        lines{end+1,1} = sprintf('%% startup (written by %s)',CORAVERSION);
        lines{end+1,1} = 'fprintf(''Startup file: %s\n'',which(mfilename))';
        lines{end+1,1} = '';
        lines{end+1,1} = '%% init tbxmanager';
        lines{end+1,1} = sprintf('addpath(genpath(''%s''))',tbxpath);
        lines{end+1,1} = 'tbxmanager restorepath';
        lines{end+1,1} = '';
        % lines{end+1,1} = '%% init MultiPrecision toolbox';
        % lines{end+1,1} = sprintf('addpath(genpath(''%s''))',mppath);
        % lines{end+1,1} = 'fprintf(''Toolbox "mp" added to the Matlab path.\n'');';
        lines{end+1,1} = '';
        lines{end+1,1} = '%% init cora';
        lines{end+1,1} = sprintf('addpath(genpath(''%s''));', CORAROOT);
        lines{end+1,1} = 'if ~isempty(which(''CORAVERSION''))';
        lines{end+1,1} = '    fprintf(''Toolbox "%s" added to the Matlab path.\n'',CORAVERSION);';
        if interactive
            % only show tip of the day for actual users, not e.g.,
            % in docker for repeatability or ci
            lines{end+1,1} = '    showCORAtipoftheday';
        end
        lines{end+1,1} = 'end';
        lines{end+1,1} = '';
        % restrict number of CPU threads
        if ~interactive
            lines{end+1,1} = '%% Restrict number of CPU threads';
            maxThreads = 4;
            lines{end+1,1} = sprintf('disp(''Limiting the number of cpu threads to %i. Change using maxNumCompThreads(N).'')', maxThreads);
            lines{end+1,1} = sprintf('maxNumCompThreads(%i);', maxThreads);
            lines{end+1,1} = '';
        end
        lines{end+1,1} = 'disp(''Done.'')';
        % finish startup function
        if isFunction
            lines{end+1,1} = '';
            lines{end+1,1} = 'end';
        end
        % end startup file
        lines{end+1,1} = '';

        % write lines back to file
        writecell(lines, startuppath, 'FileType', 'text', 'QuoteStrings', 'none');

        if interactive
            % set userpath to startup file (does not work in docker)
            startupdir = fileparts(startuppath);
            userpath(startupdir);
            addpath(startupdir);
        end

        % save paths
        addpath(tbxpath);
        addpath(currentpath);
        % savepath([userpath filesep 'pathdef.m']);
        
        % change directory back to current directory
        cd(currentpath);

        % throw error message if sedumi is missing
        if ~isempty(errortext)
            throw(CORAerror('CORA:install',errortext));
        end

        % success message
        fprintf('\nYalmip was installed successfully.\n\n')
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
    res = ~isempty(addons) && ismember(text, {addons.Name});
end

% ------------------------------ END OF CODE ------------------------------
