function writeTestSuiteSummary(testSuiteSummary)
% writeTestSuiteSummary - writes the summary of a test suite into the
%    workspace
%
% Syntax:
%    res = writeTestSuiteSummary(testSuiteSummary)
%
% Inputs:
%    writeTestSuiteSummary - struct
%
% Outputs:
%    res - logical
%
% See also: testSuiteSummary

% Authors:       Tobias Ladner
% Written:       21-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,1)
inputArgsCheck({{testSuiteSummary,'att','struct'}});

% open file
filename = [CORAROOT sprintf('/unitTests/ci/results/testSuiteSummary-%s.csv',testSuiteSummary.testSuite)];
fid = fopen(filename, 'w');

% write metadata as comments
fprintf(fid, '# date: %s\n', testSuiteSummary.date);
fprintf(fid, '# dateEnd: %s\n', testSuiteSummary.dateEnd);
fprintf(fid, '# hostsystem: %s\n', testSuiteSummary.hostsystem);
fprintf(fid, '# matlabversion: %s\n', testSuiteSummary.matlabversion);
fprintf(fid, '# coraversion: %s\n', testSuiteSummary.coraversion);
fprintf(fid, '# testSuite: %s\n', testSuiteSummary.testSuite);

% write table
CORAtable.printTable(struct2table(testSuiteSummary.results),'csv',{'s','d','d','.2f','d'},'fid',fid);

% close file
fclose(fid);

end

% ------------------------------ END OF CODE ------------------------------
