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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default values
[verbose,directory] = setDefaultValues(...
    {false,[CORAROOT filesep 'unitTests']},varargin);

% add underscore to prefix
prefix = [prefix '_'];

% list all files
files = findfiles(directory,true,prefix);

% list of currently open figures
prevFigures = get(groot, 'Children');

% number of files
nrTests = size(files,1);

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

    % extract the function name
    [~, fname] = fileparts(files(i).name);

    % save to struct
    results(i,1).fname = fname;

    try
        % print current test on command window
        if verbose
            fprintf(".. rng(%i); %s: ", testseeds(i),fname);
        end

        % set semi-random seed
        rng(testseeds(i));
        results(i,1).seed = testseeds(i);

        tic;
        % supress output of tests by usage of evalc
        if startsWith(fname,'test')
            % evaluate unit test
            [~,res] = evalc(fname);

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
        results(i,1).time = toc;

        % print result on command window
        if verbose
            if res
                fprintf('%s\n','passed');
            else
                fprintf('%s\n','failed');
            end
        end

    catch
        results(i,1).time = toc;
        res = false;
        if verbose
            fprintf('%s\n','exception');
        end
    end

    % save result to struct
    results(i,1).ok = res;

    % close opened figures within test/examples
    figures = get(groot, 'Children');
    close(figures(1:(length(figures)-length(prevFigures))))
end

% ------------------------------ END OF CODE ------------------------------
