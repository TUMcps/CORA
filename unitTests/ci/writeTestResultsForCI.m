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
%               'intlab', 'nn'
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% default format: last run
testSuite = setDefaultValues({'short'},varargin);
% check if correct identifier provided
inputArgsCheck({{testSuite,'str',{'short','long','flaky','intlab','mosek','mp','sdpt3','nn','examples','benchmarks'}}});

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
    if ~ismember(testResults.keys,testSuite)
        disp("No results for test suite with identifier '" + testSuite + "'");
        return
    else
        resultsTestSuite = testResults(testSuite);
    end
end

% read results table
results = resultsTestSuite.results; 

% result text shown in bitbucket
fileID = fopen('resultText.txt','w');
fprintf(fileID, '%d/%d tests passed.', nnz([results.ok]), length(results)); 
if sum(~[results.ok]) > 0
    fprintf(fileID, ' Failed tests. %s', strjoin({results(~[results.ok]).fname}, ', ')); 
end
fclose(fileID);

% result of test suite used as exit code
fileID = fopen('failed.txt','w');
if strcmp(testSuite,'flaky')
    fprintf(fileID, "%d", sum(~[results.ok]) > 2); 
else
    fprintf(fileID, "%d", any(~[results.ok])); 
end
fclose(fileID);

% ------------------------------ END OF CODE ------------------------------
