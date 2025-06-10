function testSuiteSummary = loadTestSuiteSummary(testSuite)
% loadTestSuiteSummary - loads the summary of a test suite into the
%    workspace
%
% Syntax:
%    testSuiteSummary = loadTestSuiteSummary(testSuite)
%
% Inputs:
%    testSuite - results of which test suite should be loaded:
%               'short' (default), 'long', 'flaky', 'mp', 'mosek', 'sdpt3',
%               'nn','examples','benchmark'
%
% Outputs:
%    testSuiteSummary - struct
%
% See also: runTestSuite

% Authors:       Tobias Ladner
% Written:       21-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(1,1)
inputArgsCheck({{testSuite,'str',{'short','long','flaky','mosek','mp','sdpt3','nn','examples','benchmarks','website','header'}}});

% check if file exist
filename = [CORAROOT sprintf('/unitTests/ci/results/testSuiteSummary-%s.csv',testSuite)];
if isempty(which(filename))
    throw(CORAerror('CORA:specialError',...
        sprintf("No data provided. Run `runTestSuite('%s')` to acquire data.",testSuite)));
end

% Initialize metadata struct
testSuiteSummary = struct();
lines = split(fileread(filename),newline);
commentlines = lines(cellfun(@(line) startsWith(line,'#'),lines));
contentlines = lines(cellfun(@(line) ~isempty(strtrim(line)) && ~startsWith(line,'#'),lines));

% Read comment lines as metadata
for i=1:numel(commentlines)
    line = commentlines{i};
    tokens = regexp(line, '#\s*(\w+):\s*(.*)', 'tokens');
    if ~isempty(tokens)
        key = tokens{1}{1};
        value = tokens{1}{2};
        testSuiteSummary.(key) = strtrim(value);
    end
end

% parse content data
data = cellfun(@(line) split(line,';'), contentlines,'UniformOutput',false);
data = [data{:}]'; % flatten
data = cellfun(@(line) strtrim(line), data, 'UniformOutput', false);

% Reassemble the full struct
testSuiteSummary.results = struct( ...
    'fname', data(2:end,1), ...
    'ok', cellfun(@(value) str2double(value)==1,data(2:end,2),'UniformOutput',false), ...
    'seed', cellfun(@(value) str2double(value),data(2:end,3),'UniformOutput',false), ...
    'time', cellfun(@(value) str2double(value),data(2:end,4),'UniformOutput',false), ...
    'checkTime', cellfun(@(value) str2double(value)==1,data(2:end,5),'UniformOutput',false) ...
);

end

% ------------------------------ END OF CODE ------------------------------
