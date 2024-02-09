function runHeaders(varargin)
% runHeaders - runs all examples from the function headers
%
% Syntax:
%    runHeaders(varargin)
%
% Inputs:
%    verbose - (optional) show workspace output or not
%    directory - (optional) change directory
%
% Outputs:
%    -

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       11-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get the original working directory
currentDirectory = pwd;

% default settings
[verbose,directory] = setDefaultValues(...
    {true,[CORAROOT filesep 'examples']},varargin);
% ... directory currently not used

% files from which we want to run the examples from their headers
files = [
    aux_addFiles('contDynamics');
    aux_addFiles('contSet');
%     aux_addFiles('converter');
    % aux_addFiles('discrDynamics');
    aux_addFiles('global');
    aux_addFiles('hybridDynamics');
    aux_addFiles('matrixSet');
];

% exclude file paths
% contSet
files = aux_exclFiles(files, ['contSet' filesep '@conPolyZono']);
% % converter
% files = aux_exclFiles(files, ['converter' filesep 'powerSystem2cora' filesep 'cases'], 'IEEE14Parameters.m');
% files = aux_exclFiles(files, ['converter' filesep 'powerSystem2cora' filesep 'cases'], 'IEEE30Parameters.m');
% global
files = aux_exclFiles(files, ['global' filesep 'thirdparty']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'combinator']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'eq_sphere_partitions']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'm2tex']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'tprod']);
files = aux_exclFiles(files, ['global' filesep 'functions'], 'Direct.m');

% check in reverse order
% files = flipud(files);

numberOfExamples = 0;
failed = {};
emptyExamples = {};

for i=1:length(files)
    file = files(i);
    file_path = [file.folder, filesep, file.name];
    filetext = fileread(file_path);

    % Extract the function name (including class if possible)
    if contains(file_path,'@')
        atPos = strfind(file_path,'@');
        % read out from after '@' and remove '.m' (-2)
        fname = file_path(atPos+1:end-2);
        % escape all backslashes
        filesepPos = strfind(fname,filesep);
        for j=1:length(filesepPos)
            % insert backslash
            fname = [fname(1:filesepPos(j)) filesep fname(filesepPos(j)+1:end)];
            % increment all indices due to added character
            filesepPos = filesepPos + 1;
        end
    else
        [~, fname] = fileparts(files(i).name);
    end

    % split text into lines
    linestext = splitlines(filetext);
    % find line with 'Example' / 'Examples'
    startline = find(contains(linestext,'Example'),1,'first');
    if isempty(startline)
        % no example
        continue
    end

    % find next part of the header which concludes the example
    endline = min([find(contains(linestext,'Reference'),1,'first'),...
        find(contains(linestext,'Other m-files required'),1,'first'),...
        find(contains(linestext,'See also'),1,'first')]);
    if isempty(endline)
        % could not determine where example ends... 'continue' for safety
        continue;
    end

    % remove percentage from beginning of each line
    exampletext = linestext(startline+1:endline-1);
    exampletext = cellfun(@(x) x(2:end),exampletext,'UniformOutput',false);
    % remove cells which are 1x0 char
    exampletext = exampletext(~cellfun(@isempty,exampletext,'UniformOutput',true));
    % if exampletext is (more or less just) "-" or "---", skip
    if any(cellfun(@(x)startsWith(x,'    -'),exampletext,'UniformOutput',true))
        emptyExamples = [emptyExamples; {fname}]; continue
    end

    % write all lines in between to text and run
    exampletext = strjoin(exampletext,'\n');

    % Supress output of example by usage of evalc
    try
        if verbose
            fprintf(['run ' fname ': ']);
        end
        % call evalc in different scope to avoid storing variables here
        aux_runSingleExample(exampletext);
        if verbose
            fprintf('%s\n','passed');
        end
    catch ME
        % save file names of failed tests
        failed = [failed; {fname}];
        if verbose
            fprintf('%s\n','failed');
        end
    end
    % close figures (if any)
    close all;

    % number of examples
    numberOfExamples = numberOfExamples + 1;

end

% display number of failed examples and file names
disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfExamples) ' examples, ' int2str(size(failed,1)) ' failed.']);
disp(strjoin(failed, ',\n'));

cd([CORAROOT filesep 'unitTests']);
save('failed.mat','failed');

% return to original working directory
cd(currentDirectory);

end


% Auxiliary functions -----------------------------------------------------

function files = aux_addFiles(path, includeSubfolders)
    if nargin < 2
        includeSubfolders = true;
    end

    subpath = '';
    if includeSubfolders
        subpath = ['**' filesep];
    end

    files = dir([CORAROOT filesep path filesep subpath '*.m']);
end

function files = aux_exclFiles(files, path, name)
    idx = contains({files.folder}', [CORAROOT filesep path]);
    if nargin == 3
        idx = idx & contains({files.name}', name);
    end
    files = files(~idx);
end

function aux_runSingleExample(exampletext)
    % run example here to have a different variable scope
    evalc(exampletext);
end

% ------------------------------ END OF CODE ------------------------------
