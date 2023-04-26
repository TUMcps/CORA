function runTestSuite(varargin)
% runTestSuite - runs a test suite by executing all functions starting with
%    the same prefix
%
% Syntax:  
%    runTestSuite()
%    runTestSuite(testSuite)
%    runTestSuite(testSuite,verbose)
%    runTestSuite(testSuite,verbose,directory)
%
% Inputs:
%    testSuite - (optional) name for test suite (case insensitive)
%                   'short' (default): prefix = 'test'
%                   'long': prefix = 'testLong'
%                   'mp': prefix = 'testMP'
%                   'mosek': prefix = 'testMOSEK'
%                   'sdpt3': prefix = 'testSDPT3'
%                   'nn': prefixes 'test_nn', 'test_neurNetContrSys', 'testnn'
%    verbose - (optional) true/false for logging results on command window
%                         default = true
%    directory - (optional) directory from which to search for tests
%                         default = cora/unitTests
%
% Outputs:
%    (text to console)
%
% Example:
%    runTestSuite();
%    runTestSuite('long',false,[CORAROOT filesep 'unitTests' filesep 'contSet' filesep 'capsule']);

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      31-August-2016
% Last update:  22-January-2021 (MW, save results of full run)
% Last revision:09-April-2023 (MW, unify all runTestSuite_* files)

%------------- BEGIN CODE --------------

% get the original working directory
currentDirectory = pwd;

% set default values
rootUnitTests = [CORAROOT filesep 'unitTests'];
[testSuite,verbose,directory] = setDefaultValues(...
    {'short',true,rootUnitTests},varargin);

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
if strcmp(directory,rootUnitTests)
    timeStart = char(datetime(now,'ConvertFrom','datenum'));
end

% find out prefix
testSuite = lower(testSuite);
switch testSuite
    case 'short'
        prefix = 'test';
    case 'long'
        prefix = 'testLong';
    case 'intlab'
        prefix = 'testINTLAB';
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
    otherwise
end

% currently, onle nn test suite has three calls
if ~iscell(prefix)
    prefix = {prefix};
end

% run main program performing the tests
results = [];
for i=1:length(prefix)
    results_i = testSuiteCore(prefix{i},verbose,directory);
    results = [results; results_i];
end

% some test suites switch the directory...
cd(directory);

% save result (only full runs)
if strcmp(directory,rootUnitTests)
    % end time
    timeEnd = char(datetime(now,'ConvertFrom','datenum'));

    % save end date and time
    data.date = [timeStart ' - ' timeEnd(13:end)];
    % save full results
    data.results = results;

    % save matlab version
    data.matlabversion = version;

    % save operating system
    data.hostsystem = 'Unknown';
    if ismac
        data.hostsystem = 'Mac';
    elseif isunix
        data.hostsystem = 'Unix';
    elseif ispc
        data.hostsystem = 'Windows';
    end

    % use map to save results (we use containers.Map for legacy support
    % over dictionary which only exists since R2022b)
    newResults = containers.Map(testSuite,data);

    try
        % load data from .mat file
        load('unitTestsStatus.mat','testResults');
        % append result to map (if the key has already been in use, the new
        % results overwrite the old results)
        testResults = [testResults; newResults];
    catch
        % no previous results saved -> save directly to .mat file
        testResults = newResults;
    end
    % save results
    save('unitTestsStatus.mat','testResults');

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

%------------- END OF CODE --------------