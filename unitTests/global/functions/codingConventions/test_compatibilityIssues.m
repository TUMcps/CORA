function res = test_compatibilityIssues()
% test_compatibilityIssues - tests issues found by the MATLAB Code Analyzer
%
% Syntax:
%    res = test_compatibilityIssues()
%
% Inputs:
%    -
%
% Outputs:
%    res - whether all files are valid
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: codeCompatibilityReport, codeIssues

% Authors:       Tobias Ladner
% Written:       05-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% run MATLAB Code Analyzer and extract issues
fprintf('Analyzing all files in %s..',CORAROOT)
issues = codeIssues(CORAROOT);
issues = issues.Issues;
fprintf('\n\n')

% reduce severity for certain problems as the analyzer is doesn't care
% about 'isMATLABReleaseOlderThan'
issues{contains(issues.Description,'dual-simplex-legacy'),2} = matlab.codeanalysis.IssueSeverity.warning;

% categorize issues
errorIssues = issues(issues.Severity == matlab.codeanalysis.IssueSeverity.error,:);
warningIssues = issues(issues.Severity == matlab.codeanalysis.IssueSeverity.warning,:);
infoIssues = issues(issues.Severity == matlab.codeanalysis.IssueSeverity.info,:);

% print fatal issues
if height(errorIssues) > 0
    disp('Urgent issues:')
    CORAtable.printTable(errorIssues,'single');
end
fprintf('%i urgent issues found. %i warnings. %i infos. Open <a href="matlab:codeCompatibilityReport">Code Compatibility Analyzer</a> to see all issues.\n', height(errorIssues),height(warningIssues),height(infoIssues))

% get results
res = isempty(errorIssues);

end

% ------------------------------ END OF CODE ------------------------------
