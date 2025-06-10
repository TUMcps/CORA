function runTestSuite(varargin)
% runTestSuite - runs a test suite by executing all functions starting with
%    the same prefix
%
% Syntax:
%    runTestSuite()
%    runTestSuite(testSuite)
%    runTestSuite(testSuite,verbose)
%    runTestSuite(testSuite,verbose,directory,checkTime)
%
% Inputs:
%    testSuite - (optional) name for test suite (case insensitive)
%                   'short' (default): prefix = 'test'
%                   'long': prefix = 'testLong'
%                   'nn': prefixes 'test_nn', 'test_neurNetContrSys', 'testnn'
%                   'examples': prefix 'example'
%                   'benchmark': prefix 'benchmark'
%                   'website': prefix 'website'
%                   'testHeader': examples in docstring headers
%                   'mp': prefix = 'testMP'
%                   'mosek': prefix = 'testMOSEK'
%                   'sdpt3': prefix = 'testSDPT3'
%                 adding ':failed', only tests that failed in the last run
%                 are executed
%    verbose - (optional) true/false for logging results on command window
%                         default = true
%    directory - (optional) directory from which to search for tests
%                         default = cora/unitTests
%    checkTime - (optional) logical, compare execution time with saved run
%
% Outputs:
%    (text to console)
%
% Example:
%    runTestSuite();
%    runTestSuite('long',false,[CORAROOT filesep 'unitTests' filesep 'contSet' filesep 'capsule']);

% Authors:       Matthias Althoff, Mark Wetzlinger, Tobias Ladner
% Written:       31-August-2016
% Last update:   22-January-2021 (MW, save results of full run)
%                16-June-2023 (MW, add ':failed')
%                21-May-2025 (TL, print results instead of .mat)
% Last revision: 09-April-2023 (MW, unify all runTestSuite_* files)

% ------------------------------ BEGIN CODE -------------------------------

% get the original working directory
currentDirectory = pwd;

% set default values
[testSuite,verbose,directory,checkTime] = setDefaultValues(...
    {'short',true,[],true},varargin);

% set directory
rootUnitTests = [CORAROOT filesep 'unitTests'];
rootExamples = [CORAROOT filesep 'examples'];
if isempty(directory) 
    switch testSuite
        case {'examples','benchmarks','website'}
            directory = rootExamples;
        otherwise
            directory = rootUnitTests;
    end
end

% correct separator for directory
temp = strrep(directory,'\',filesep);
temp = strrep(temp,'/',filesep);
directory = temp;
% check if directory actually exists
if ~isfolder(directory)
    throw(CORAerror('CORA:specialError',...
        "Input given for directory is not a directory!"));
end

% save start date and time (only full runs)
isFullRun = strcmp(directory,rootUnitTests) || strcmp(directory,rootExamples);
if isFullRun
    timeStart = datetime('now','Format','dd-MMM-yyyy HH:mm:ss');
end

% robustify
testSuite = lower(testSuite);

% check whether only failed tests are wanted
if endsWith(testSuite,':failed')
    % remove ':failed' from char array
    testSuite = testSuite(1:end-7);
    aux_runFailedTests(testSuite);
    return
end

% find out prefix
switch testSuite
    case 'short'
        prefix = 'test';
    case 'long'
        prefix = 'testLong';
    case 'flaky'
        prefix = 'testFlaky';
    case 'mosek'
        prefix = 'testMOSEK';
    case 'mp'
        prefix = 'testMP';
    case 'nn'
        % tests run with default CORA installation
        prefix{1} = 'test_nn';
        prefix{2} = 'test_neurNetContrSys';
        % tests require Deep Learning Toolbox + ONNX Converter
        prefix{3} = 'testnn';
    case 'sdpt3'
        prefix = 'testSDPT3';
    case 'examples'
        prefix = 'example';
    case 'benchmarks'
        prefix = 'benchmark';
    case 'website'
        prefix = 'website';
    case 'header'
        prefix = 'testHeader';
    otherwise
        throw(CORAerror('CORA:specialError','Unknown test suite.'))
end
% convert to cell. Currently, only nn test suite has three calls
if ~iscell(prefix)
    prefix = {prefix};
end

% load previous results (or init empty)
oldResults = [];
if checkTime
    try
        oldTestSuiteSummary = loadTestSuiteSummary(testSuite);
        oldResults = oldTestSuiteSummary.results;
        disp('Previous results loaded successfully.')
    catch ME
        fprintf('! No previous results loaded (%s).\n',ME.message)
    end
else
    disp('Not checking previous results.')
end

% run main program performing the tests
results = [];
for i=1:length(prefix)
    results_i = testSuiteCore(prefix{i},verbose,directory,oldResults);
    results = [results; results_i];
end

% some test suites switch the directory...
cd(currentDirectory);

% save result (only full runs)
if isFullRun
    % init results struct
    testSuiteSummary = struct;

    % save end date and time
    testSuiteSummary.date = timeStart;
    testSuiteSummary.dateEnd = datetime('now','Format','dd-MMM-yyyy HH:mm:ss');

    % save operating system
    testSuiteSummary.hostsystem = 'Unknown';
    if ismac
        testSuiteSummary.hostsystem = 'Mac';
    elseif isunix
        testSuiteSummary.hostsystem = 'Unix';
    elseif ispc
        testSuiteSummary.hostsystem = 'Windows';
    end

    % save matlab version
    testSuiteSummary.matlabversion = version;

    % save CORA version
    testSuiteSummary.coraversion = CORAVERSION;

    % save test suite
    testSuiteSummary.testSuite = testSuite;
    
    % save full results
    testSuiteSummary.results = results;
    
    % save test suite summary
    writeTestSuiteSummary(testSuiteSummary)

    % print result
    printTestOverview(testSuite);
    
else
    % display number of failed tests + file names
    disp('----------------------------------------------------------------------------');
    fprintf("run %i tests, %i failed.\n",length(results),nnz(~[results.ok]));
    disp(strjoin("  " + {results(~[results.ok]).fname}, ',\n'));
    
end

% return to original working directory
cd(currentDirectory);

end


% Auxiliary functions -----------------------------------------------------

function aux_runFailedTests(testSuite)
% executes all unit tests of a given test suite that failed in the last
% full run of that test suite

    % load data from .mat file
    try
        load('unitTestsStatus.mat','testResults');
    catch
        % no previous results saved
        disp("No saved results available! " + ...
            "Run 'runTestSuite('" + testSuite + "')'; to generate.");
        return
    end

    % read out failed tests
    res = testResults(testSuite).results;

    % list of failed tests
    allTests = {res.fname}';
    failedTests = allTests(~logical([res.ok]'));

    % execute tests
    for i=1:length(failedTests)
        % read out name
        fname = failedTests{i};

        try
            fprintf(['run ' fname ': ']);
            [~,res] = evalc(fname);
            if res
                fprintf('%s\n','passed');
            else
                fprintf('%s\n','failed');
            end
        catch
            fprintf('%s\n','failed');
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
