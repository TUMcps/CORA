function printTestOverview(varargin)
% printTestOverview - print information about the current state
%    of unit testing to the console
%
% Syntax:  
%    printTestOverview()
%    printTestOverview(format)
%    printTestOverview(format,')
%
% Inputs:
%    format - (optional) results of which test suite should be displayed:
%               'short' (default), 'long', 'mp', 'mosek', 'sdpt3',
%               'intlab', 'nn'
%    type - name-value pairs
%               'LongestTests', <longestTests> number of longest tests
%                                              (default: 5)
%
% Outputs:
%    -
%
% Example: 
%    printTestOverview('long');
%    printTestOverview('long','LongestTests',10);

% Author:       Mark Wetzlinger
% Written:      22-January-2021
% Last update:  09-April-2023 (MW, remove 'classes', integrate other test suites)
%               28-April-2023 (MW, add name-value pairs)
% Last revision:---

%------------- BEGIN CODE --------------

if nargin > 1 && mod(nargin-1,2) ~= 0
    throw(CORAerror('CORA:oddNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(2:end);
end

% default format
format = setDefaultValues({'short'},varargin);
% more leniency in writing of identifiers
format = lower(format);
% check if correct identifier provided
inputArgsCheck({{format,'str',{'short','long','intlab','mosek','mp','sdpt3','nn'}}});

% check list of name-value pairs
checkNameValuePairs(NVpairs,{'LongestTests'});
% read name-value pairs
[NVpairs,longestTests] = readNameValuePair(NVpairs,'LongestTests','isscalar',5);


% load data
unitTestsFile = [CORAROOT filesep 'unitTests' filesep 'unitTestsStatus.mat'];
if ~isfile(unitTestsFile)
    throw(CORAerror('CORA:specialError',...
        "No data provided. Run '''runTestSuite''' to acquire data."));
else
    % latest test results are stored in the variable 'testResults' in
    % unitTestsStatus.mat
    load(unitTestsFile,'testResults');
    % read out data from map
    if ~ismember(testResults.keys,format)
        disp("No results for test suite with identifier '" + format + "'");
        return
    else
        resultsTestSuite = testResults(format);
    end
end

% save results to shorter variable
results = resultsTestSuite.results;

% cap value for number of longest tests
if longestTests > length(results)
    longestTests = length(results);
end

% read data from results
nrTotalTests = length(results);
nrFailedTests = nnz(~[results.ok]);
nameFailedTests = {results(~[results.ok]).fname};
seedFailedTests = [results(~[results.ok]).seed];

% print header
fprintf('-*---------------------------------*-\n');
fprintf('--- Results of last full test run ---\n\n');

% print which test suite
fprintf(['  test suite: ',format,'\n']);

% print machine and Matlab version
fprintf(['  system:     ',resultsTestSuite.hostsystem,'\n']);
fprintf(['  Matlab:     ',resultsTestSuite.matlabversion,'\n']);

% print date and time
fprintf(['  date/time:  ',resultsTestSuite.date,'\n']);

% print test results
fprintf('  tests:\n');
fprintf(['  .. total:   ',num2str(nrTotalTests),'\n']);
fprintf(['  .. failed:  ',num2str(nrFailedTests),'\n']);

% print names of all failed tests
if nrFailedTests > 0
    fprintf("  .. .. %s (seed=%i)\n", nameFailedTests{1}, seedFailedTests(1));
    for i=2:nrFailedTests
    fprintf("        %s (seed=%i)\n", nameFailedTests{i}, seedFailedTests(i));
    end
end

% read out longest tests
if longestTests > 0
    [compTime,idxMax] = maxk([results.time],longestTests);
    fprintf("  .. longest:\n");
    for i=1:longestTests
    fprintf("  .. .. %s (%.2fs)\n",results(idxMax(i)).fname,compTime(i));
    end
end

% print footer
fprintf('-*---------------------------------*-\n');

%------------- END OF CODE --------------