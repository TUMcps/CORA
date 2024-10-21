function writeTestResultsForCI(varargin)
% writeTestResultsForCI - prepares test result for gitlab CI
%
% Syntax:
%    writeTestResultsForCI()
%    writeTestResultsForCI(testSuite)
%
% Inputs:
%    testSuite - (optional) results of which test suite should be displayed:
%               'short' (default), 'long', 'mp', 'mosek', 'sdpt3',
%               'intlab', 'nn', 'examples','benchmarks','header'
%
% Outputs:
%    -
%
% Example: 
%    writeTestResultsForCI('long');
%
% See also: printTestOverview

% Authors:       Tobias Ladner
% Written:       22-April-2023
% Last update:   10-May-2024 (TL, removed resultText.txt for bitbucket)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume unsuccessful test suite
failed = true;

% default format: last run
testSuite = setDefaultValues({'short'},varargin);
% check if correct identifier provided
inputArgsCheck({{testSuite,'str',{'short','long','flaky','intlab','mosek','mp','sdpt3','nn','examples','benchmarks','header'}}});

% load data
unitTestsFile = [CORAROOT filesep 'unitTests' filesep 'unitTestsStatus.mat'];
if ~isfile(unitTestsFile)
    disp("No data provided. Run '''runTestSuite''' to acquire data.");
else
    % latest test results are stored in the variable 'testResults' in
    % unitTestsStatus.mat
    load(unitTestsFile,'testResults');
    % read out data from map
    if ~ismember(testResults.keys,testSuite)
        disp("No results for test suite with identifier '" + testSuite + "'");
    else
        % read results table
        resultsTestSuite = testResults(testSuite);
        results = resultsTestSuite.results; 
        if strcmp(testSuite,'flaky')
            failed = sum(~[results.ok]) > 2; 
        else
            failed = any(~[results.ok]); 
        end
    end
end

% result of test suite used as exit code
fileID = fopen('failed.txt','w');
% EXIT_CODE = 0 means successful test suite
fprintf(fileID, "%d", failed); 
fclose(fileID);

end

% ------------------------------ END OF CODE ------------------------------
