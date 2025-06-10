function writeTestResultsForCI(varargin)
% writeTestResultsForCI - prepares test result for gitlab CI
%    in particular, writes a file EXIT_CODE.txt containing the exit code
%    for the job. See .gitlab-ci.yml for full list of exit codes.
%
% Syntax:
%    writeTestResultsForCI()
%    writeTestResultsForCI(testSuite)
%
% Inputs:
%    testSuite - (optional) results of which test suite should be displayed:
%               'short' (default), 'long', 'mp', 'mosek', 'sdpt3',
%               'nn', 'examples','benchmarks','header'
%
% Outputs:
%    -
%
% See also: .gitlab-ci.yml, printTestOverview

% Authors:       Tobias Ladner
% Written:       22-April-2023
% Last update:   10-May-2024 (TL, removed resultText.txt for bitbucket)
%                11-April-2025 (TL, proper EXIT_CODES)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume EXIT_SCORE=5 (unknown error)
EXIT_CODE = 5;

% load test suite
testSuite = setDefaultValues({'short'},varargin);
testSuiteSummary = loadTestSuiteSummary(testSuite);

% check if test results are recent (<= 5min)
if minutes(datetime('now')-datetime(testSuiteSummary.dateEnd)) <= 5
    % read results table
    results = testSuiteSummary.results;
    % jobs succeeds with EXIT_CODE=0, and fails for EXIT_CODE=1
    if strcmp(testSuite,'flaky')
        EXIT_CODE = sum(~[results.ok]) > 2; 
    else
        EXIT_CODE = any(~[results.ok]); 
    end
end

% result of test suite used as exit code
fileID = fopen('EXIT_CODE.txt','w');
fprintf(fileID, "%d", EXIT_CODE); 
fclose(fileID);

end

% ------------------------------ END OF CODE ------------------------------
