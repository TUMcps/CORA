function results = testSuiteCore(prefix,varargin)
% testSuiteCore - runs functions starting with a certain prefix contained 
%    in the directory and recursively searches all subfolders for same prefix
%
% Syntax:
%    results = testSuiteCore(prefix)
%    results = testSuiteCore(prefix,verbose,directory,oldResults)
%
% Inputs:
%    prefix - prefix of function names to be tested
%    verbose - (optional) show workspace output or not
%    directory - (optional) change directory
%    oldResults - (optional) struct containing old results
%
% Outputs:
%    results - struct with results of test suite, fields
%              'fname': function name of test
%              'ok': true/false whether test successful
%              'seed': randomized seed for test
%              'time': computation time for test
%
% Example:
%    -

% Authors:       Dmitry Grebenyuk, Matthias Althoff, Tobias Ladner
% Written:       31-August-2016
% Last update:   22-April-2024 (TL, better formatting)
%                25-September-2024 (TL, more expressive error messages)
%                19-May-2025 (TL, added option to compare time w/ prev. runs)
%                23-May-2025 (TL, pre-init results struct, deal w/ parallel pool)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[verbose,directory,oldResults] = setDefaultValues(...
    {false,[CORAROOT filesep 'unitTests'],[]},varargin);

% add underscore to prefix
prefix = [prefix '_'];

% file extension
fileext = 'm';
if strcmp(prefix,'website_')
    fileext = 'mlx';
end

% list all files
files = findfiles(directory,true,prefix,fileext);
% filter s.t. test_codingConventions is always on top (if present)
idxCC = arrayfun(@(file) strcmp(file.name,'test_codingConventions.m'), files);
files = [files(idxCC); files(~idxCC)];

% number of files
nrTests = size(files,1);
maxFileName = max(arrayfun(@(file) length(file.name)-2, files));

% check for uniqueness of filename
if numel(unique({files.name})) ~= nrTests
    aux_handleNonUniqueFileNames(files);
end

% generate random seed for entire test suite from current time
% (required for CI to run different seeds every time)
rng(convertTo(datetime(),'posixtime'))

% init struct
results = struct( ...
    'fname',arrayfun(@(file) file.name(1:(end-numel(fileext)-1)),files,'UniformOutput',false), ...
    'ok',num2cell(nan(nrTests,1)), ...
    'seed',num2cell(randi(2^32,nrTests,1)), ...
    'time',num2cell(nan(nrTests,1)), ...
    'checkTime',num2cell(true(nrTests,1)) ...
);

% check surrounding ---
% list of currently open figures
prevFigures = get(groot, 'Children');
% check if parallel computing toolbox is installed
hasPCT = ~isempty(which('gcp'));

% loop over all test files
fprintf("Running %g tests:\n", nrTests)
for i=1:nrTests
    % read out properites
    fname = results(i,1).fname;
    seed = results(i,1).seed;

    try
        % print current test on command window
        if verbose
            fprintf(['.. rng(%10i); %-' num2str(maxFileName) 's : '],seed,fname);
        end

        % set semi-random seed
        rng(seed);

        testTime = tic;
        % suppress output of tests by usage of evalc 
        % except for the header test, as it has its own output
        if startsWith(fname, 'testHeader')
            res = eval(fname);
        elseif startsWith(fname,'test')
            % evaluate unit test
            [~,res] = evalc(fname);
        elseif startsWith(fname,'website')
            % live scripts cannot be executed directly with evalc and leak
            % variables outside of their scope. 
            % -> put them into a function handle to suppress outputs
            % and do not overwrite any variables in current workspace
            aux_mywebsitefun(fname);            
            res = true;
        else
            % no output required for examples/benchmarks/...
            try
                evalc(fname);
            catch ME
                if strcmp(ME.identifier, 'CORA:noSuitableSolver') ...
                   % missing solver is ok
                else
                    rethrow(ME)                    
                end
            end
            res = true;
        end

        % save time
        results(i,1).time = toc(testTime);

        % print result on command window
        if verbose
            if res
                % check if it took much longer than previous runs
                msg = '';
                if ~isempty(oldResults)
                    % check if current test is present
                    oldTestResult = oldResults(strcmp({oldResults.fname},fname));
                    if isempty(oldTestResult)
                        % no time saved
                        msg = '! No previous time saved.';
                    elseif oldTestResult.checkTime
                        % read out time
                        oldTime = oldTestResult.time;
                        newTime = results(i,1).time;
                        % check if time difference is reasonable
                        if (2*oldTime+3 < newTime) || ((oldTime-1)/3 > newTime)
                            % test took too long/short, worth double-checking...
                            msg = sprintf( ...
                                '! Test passed but took unexpectedly long/short: %.2fs -> %.2fs (%+.0f%%).', ...
                                oldTime,newTime,((newTime/oldTime)-1)*100 ...
                            );
                        end
                    else
                        % also don't check time in the future
                        results(i,1).checkTime = false;
                    end
                end
                if ~isempty(msg)
                    % add formatting
                    msg = sprintf(' : %s',msg);
                end
                fprintf('%s  (%.2fs)%s\n','passed',results(i,1).time,msg);
            else
                fprintf('%s  (%.2fs) : returned ''false''\n','failed',results(i,1).time);
            end
        end

    catch ME
        results(i,1).time = toc(testTime);
        res = false;
        if verbose
            % read error message and line number
            idxThisFile = find(arrayfun(@(entry) contains(entry.name,mfilename), ME.stack),1,'last')';
            idxTest = max([1,idxThisFile-1]);
            entryTestFile = ME.stack(idxTest);
            msg = strrep(ME.message,newline,' | '); % single line

            % read file where error is thrown for additional information
            addInfo = "";
            % shift to file where error is thrown
            idxErrorFile = 1;
            while ismember(ME.stack(idxErrorFile).name, {'assertLoop'})
                % go to 'real' file where error is thrown
                idxErrorFile = idxErrorFile + 1;
            end
            if idxErrorFile ~= idxTest
                entryErrorFile = ME.stack(idxErrorFile);
                filepathErrorFile = entryErrorFile.file;
                if contains(filepathErrorFile,CORAROOT)
                    filepathErrorFile = filepathErrorFile(numel(CORAROOT)+2:end);
                end
                addInfo = sprintf(' (thrown in %s, line %i)',filepathErrorFile,entryErrorFile.line);
            end

            % display results
            fprintf('%s  (%.2fs) : in line %i: ''%s''%s\n','failed',results(i,1).time,entryTestFile.line,msg,addInfo);
        end
    end

    % save result to struct
    results(i,1).ok = res;

    % close opened figures within test/examples
    figures = get(groot, 'Children');
    close(figures(1:(length(figures)-length(prevFigures))))

    % shut down parallel pool to not interfer with other tests
    if hasPCT
        delete(gcp("nocreate"))
    end
end

end


% Auxiliary functions -----------------------------------------------------

function aux_handleNonUniqueFileNames(files)

% Find unique strings and index mapping
[uniqueStrs, ~, idx] = unique({files.name});

% Count occurrences of each unique string
counts = histcounts(idx, 1:(max(idx)+1));

% Find duplicates (indices in 'strs' that are duplicates)
dupidx = find(ismember(idx, find(counts > 1)))';

% display duplicates
disp('ERROR: Found non-unique filenames!')
for i=dupidx
    fprintf('- %s/%s\n',files(i).folder,files(i).name);
end

% throw error as matlab will always run the last test on its path, which is
% undesired behavior
throw(CORAerror('CORA:notSupported', ...
    ['Unable to run test suite where some tests have the same name.\n' ...
    'This would result in always runnig the last test on the path, leaving all former tests untested.\n' ...
    'Please either merge the two tests or give them unique names to better state their purpose.']))

end

function aux_mywebsitefun(fname)
    % this auxiliary function prevents variables within live scripts to
    % leak out into the workspace
    evalc('run(fname)');
end

% ------------------------------ END OF CODE ------------------------------
