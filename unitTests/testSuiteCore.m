function results = testSuiteCore(prefix,varargin)
% testSuiteCore - runs functions starting with a certain prefix contained 
%    in the directory and recursively searches all subfolders for same prefix
%
% Syntax:
%    results = testSuiteCore(prefix)
%    results = testSuiteCore(prefix,verbose)
%    results = testSuiteCore(prefix,verbose,directory)
%
% Inputs:
%    prefix - prefix of function names to be tested
%    verbose - (optional) show workspace output or not
%    directory - (optional) change directory
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

% Authors:       Dmitry Grebenyuk, Matthias Althoff
% Written:       31-August-2016
% Last update:   22-April-2024 (TL, better formatting)
%                25-September-2024 (TL, more expressive error messages)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[verbose,directory] = setDefaultValues(...
    {false,[CORAROOT filesep 'unitTests']},varargin);

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

% list of currently open figures
prevFigures = get(groot, 'Children');

% number of files
nrTests = size(files,1);
maxFileName = max(arrayfun(@(file) length(file.name)-2, files));

% generate random testSuite seed from current time
% (required for CI to run different seeds every time)
rng(convertTo(datetime(),'posixtime'))

% generate test seeds
testseeds = randi(2^32,nrTests,1);

% init empty struct
results = struct('fname',{},'ok',{},'seed',{},'time',{});

% loop over all test files
fprintf("Running %g tests:\n", nrTests)
for i=1:nrTests

    % extract the function name and path
    [~, fname] = fileparts(files(i).name);
    fpath = [files(i).folder filesep files(i).name];    

    % save to struct
    results(i,1).fname = fname;

    try
        % print current test on command window
        if verbose
            fprintf(['.. rng(%10i); %-' num2str(maxFileName) 's : '], ...
                testseeds(i),fname);
        end

        % set semi-random seed
        rng(testseeds(i));
        results(i,1).seed = testseeds(i);

        testTime = tic;
        % suppress output of tests by usage of evalc 
        % except for the header test, as it has its own output
        if startsWith(fname, 'testHeader')
            res = eval(fname);
        elseif startsWith(fname,'test')
            % evaluate unit test
            [~,res] = evalc(fname);
        elseif startsWith(fname,'website')
            % live scripts cannot be executed with evalc
            % -> put them into a function handle to suppress outputs
            evalc('run(fpath);');
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
                fprintf('%s  (%.2fs)\n','passed',results(i,1).time);
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
end

% ------------------------------ END OF CODE ------------------------------
