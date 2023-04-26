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

% Author:       Dmitry Grebenyuk, Matthias Althoff
% Written:      31-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values
[verbose,directory] = setDefaultValues(...
    {false,[CORAROOT filesep 'unitTests']},varargin);

% add underscore to prefix
prefix = [prefix '_'];

% list all files
files = findfiles(directory,true,prefix);

% init struct
results = struct('fname',[],'ok',[],'seed',[],'time',[]);

% number of files
nrTests = size(files,1);
% generate test seeds
testseeds = randi(2^32,nrTests,1);

% loop over all test files
for i=1:nrTests

    % extract the function name
    [~, fname] = fileparts(files(i).name);

    % save to struct
    results(i,1).fname = fname;

    try
        % print current test on command window
        if verbose
            fprintf(['run ' fname ': ']);
        end

        % set semi-random seed
        rng(testseeds(i));
        results(i,1).seed = testseeds(i);

        tic;
        % supress output of tests by usage of evalc
        [~,res] = evalc(fname);

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
            fprintf('%s\n','failed');
        end
    end

    % save result to struct
    results(i,1).ok = res;
end

%------------- END OF CODE --------------