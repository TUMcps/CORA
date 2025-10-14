function res = runHeaders(varargin)
% runHeaders - runs all examples from the function headers
%
% Syntax:
%    runHeaders(varargin)
%
% Inputs:
%    verbose - (optional) show workspace output or not
%
% Outputs:
%    res - if all headers ran successfully

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       11-April-2023
% Last update:   04-March-2025 (TL, all folders are tested now)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get the original working directory
currentDirectory = pwd;

% default settings
[verbose] = setDefaultValues({true},varargin);

% files from which we want to run the examples from their headers
files = [
    aux_addFiles('app');
    aux_addFiles('contDynamics');
    aux_addFiles('contSet');
    aux_addFiles('converter');
    aux_addFiles('discrDynamics');
    aux_addFiles('global');
    aux_addFiles('examples');
    aux_addFiles('hybridDynamics');
    aux_addFiles('matrixSet');
    aux_addFiles('nn');
    aux_addFiles('specification');
    aux_addFiles('unitTests');
];

% exclude:
% - thirdparty files
files = aux_exclFiles(files, ['global' filesep 'thirdparty']);

% check in reverse order
% files = flipud(files);

numberOfExamples = 0;
failed = {};
emptyExamples = {};

fprintf("Testing up to %i headers:\n", length(files))
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
            fname = [fname(1:filesepPos(j)-1) '/' fname(filesepPos(j)+1:end)];
        end
    else
        [~, fname] = fileparts(files(i).name);
    end

    % split text into lines
    linestext = splitlines(filetext);
    % find line with 'Example' / 'Examples'
    startline = find(contains(linestext,['%' ' Example']),1,'first');
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

    seed = randi(2^32);

    % Suppress output of example by usage of evalc
    try
        if verbose
            fprintf('.. rng(%10i); header of %-70s : ', seed, fname);
        end
        % call evalc in different scope to avoid storing variables here
        aux_runSingleExample(exampletext, seed);
        if verbose
            fprintf('%s\n','passed');
        end
    catch ME
        if contains(file_path,[filesep 'private' filesep])
            % ok as private functions cannot be called from outside
            if verbose
                fprintf('%s\n','maybe');
            end
        else
            % save file names of failed tests
            failed = [failed; {fname}];
            if verbose
                fprintf('%s\n','failed');
            end
        end
    end
    % close figures (if any)
    close all;

    % number of examples
    numberOfExamples = numberOfExamples + 1;

end

numFailed = size(failed,1);

% display number of failed examples and file names
disp('----------------------------------------------------------------------------');
fprintf('run %i examples in headers, %i failed.\n', numberOfExamples, numFailed);
if numFailed > 0
    fprintf('- header of ')
    disp(strjoin(failed, '\n- header of '));
end

cd([CORAROOT filesep 'unitTests']);
save('failed.mat','failed');

% return to original working directory
cd(currentDirectory);

% check if all succeeded
res = isempty(failed);

end


% Auxiliary functions -----------------------------------------------------

function files = aux_addFiles(path, includeSubfolders)
    % add files based on pattern
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
    % exclude files based on pattern
    idx = contains({files.folder}', [CORAROOT filesep path]);
    if nargin == 3
        idx = idx & contains({files.name}', name);
    end
    files = files(~idx);
end

function aux_runSingleExample(exampletext,seed)
    % run example here to have a different variable scope
    rng(seed)
    evalc(exampletext);
end

% ------------------------------ END OF CODE ------------------------------
