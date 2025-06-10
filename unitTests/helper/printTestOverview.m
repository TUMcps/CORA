function printTestOverview(varargin)
% printTestOverview - print information about the current state
%    of unit testing to the console
%
% Syntax:
%    printTestOverview()
%    printTestOverview(testSuite)
%    printTestOverview(format,'LongestTests',5)
%
% Inputs:
%    testSuite - (optional) results of which test suite should be displayed:
%               'short' (default), 'long', 'flaky', 'mp', 'mosek', 'sdpt3',
%               'nn'
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

% Authors:       Mark Wetzlinger
% Written:       22-January-2021
% Last update:   09-April-2023 (MW, remove 'classes', integrate other test suites)
%                28-April-2023 (MW, add name-value pairs)
%                23-November-2023 (TL, failed tests show call syntax)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin > 1 && mod(nargin-1,2) ~= 0
    throw(CORAerror('CORA:oddNumberInputArgs'));
else
    % read input arguments
    NVpairs = varargin(2:end);
end

% load test suite
testSuite = setDefaultValues({'short'},varargin);
testSuiteSummary = loadTestSuiteSummary(testSuite);

% check list of name-value pairs
checkNameValuePairs(NVpairs,{'LongestTests'});
% read name-value pairs
[NVpairs,longestTests] = readNameValuePair(NVpairs,'LongestTests','isscalar',5);

% save results to shorten variable
results = testSuiteSummary.results;

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
fprintf(['  Test suite: ',testSuite,'\n']);

% print CORA version, Matlab version and machine
if isfield(testSuiteSummary,'coraversion')
    fprintf(['  Version:    ',testSuiteSummary.coraversion,'\n']);
end
fprintf(['  Matlab:     ',testSuiteSummary.matlabversion,'\n']);
fprintf(['  System:     ',testSuiteSummary.hostsystem,'\n']);

% print date and time
fprintf(['  Date:       ',char(testSuiteSummary.date),' - ',char(testSuiteSummary.dateEnd),'\n']);

% print test results
fprintf('  Tests:\n');
fprintf(['  .. total:   ',num2str(nrTotalTests),'\n']);
fprintf(['  .. failed:  ',num2str(nrFailedTests),'\n']);

% print names of all failed tests
if nrFailedTests > 0
    fprintf("  .. .. rng(%i); %s\n", seedFailedTests(1),nameFailedTests{1});
    for i=2:nrFailedTests
    fprintf("        rng(%i); %s\n", seedFailedTests(i),nameFailedTests{i});
    end
end

% read out longest tests
if longestTests > 0
    [compTime,idxMax] = maxk([results.time],longestTests);
    fprintf("  .. longest:\n");
    for i=1:min(longestTests,numel(idxMax))
    fprintf("  .. .. %s (%.2fs)\n",results(idxMax(i)).fname,compTime(i));
    end
end

% print footer
fprintf('-*---------------------------------*-\n');

% ------------------------------ END OF CODE ------------------------------
