function res = unitTestsOverview(varargin)
% unitTestsOverview - go through cora repository to find for which functions
%    do unit tests exist
%
% Syntax:  
%    lastFullTestRun
%
% Inputs:
%    testrun - (struct) result of full test suite
%
% Outputs:
%    res - execution successful
%
% Example: 
%    -

% Author:       Mark Wetzlinger
% Written:      20-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default output argument
res = false;

% default input argument
newtestrun = false;
if nargin == 1
    if isstruct(varargin{1})
        testrun = varargin{1};
        newtestrun = true;
    end
end

% get the original working directory
currentDirectory = pwd;
% switch to directory of unit tests
directory = [coraroot filesep 'unitTests'];
cd(directory);

% check if file available (as input argument or saved .mat file)
if ~isfile('unitTestsStatus.mat') && ~newtestrun
    error("No data provided. Run '''runTestSuite''' to acquire data.");
elseif newtestrun
    % new data
    overview.lastRun = testrun;
else
    % load results from last full test run
    data = load('unitTestsStatus.mat');
    overview = data.overview;
end


% list main content folders
foldernames = {'contDynamics';'contSet';'discrDynamics';'hybridDynamics';'matrixSet'};

% write all funcs in one list and find matching unit tests
allfuncs = {};
alltests = {};
for i=1:length(foldernames)
    % list functions and unit tests
    tempfiles = listFolderContent([coraroot filesep foldernames{i}]);
    temptests = listFolderContent([coraroot filesep 'unitTests' filesep foldernames{i}]);
    % collapse lists into one list
    tempfiles = functionList(tempfiles);
    allfuncs = [allfuncs; tempfiles];
    temptests = functionList(temptests);
    alltests = [alltests; temptests];
end

corarootlength = length(coraroot);
% crop coraroot from all function names
for i=1:length(allfuncs)
    funcname = extractBetween(allfuncs{i},corarootlength+2,length(allfuncs{i}));
    allfuncs{i} = funcname{1};
end

% extract function names of each test
for i=1:length(alltests)
    idxfilesep = strfind(alltests{i},filesep);
    funcname = extractBetween(alltests{i},idxfilesep(end)+1,length(alltests{i}));
    alltests{i} = funcname{1};
end

overview.fullList.func = allfuncs;
% find unit test for each function
for i=1:length(allfuncs)
    overview.fullList.test{i,1} = findUnitTest(allfuncs{i},alltests);
    % check if run was successful in last run
    if strcmp(overview.fullList.test{i},'---')
        overview.fullList.status(i,1) = -1; % no information
    else
        if any(ismember(overview.lastRun.nameFailedTests,overview.fullList.test{i}))
            overview.fullList.status(i,1) = 0; % at least one test failed
        else
            overview.fullList.status(i,1) = 1; % all tests successful
        end
        % delete test(s) from alltests
        for j=1:length(overview.fullList.test{i})
            alltests = alltests(~strcmp(alltests,overview.fullList.test{i}{j}));
        end
        overview.fullList.test{i} = strjoin(overview.fullList.test{i},'\n');
    end
end

% unassigned test: either duplicates (e.g. reach_02) or name not found in func
overview.unassignedList.test = alltests;
if ~isempty(overview.unassignedList.test)
overview.unassignedList.status = ...
    ~ismember(overview.unassignedList.test,overview.lastRun.nameFailedTests);
end
% issue: if method from superclass called, function test in unassignedList
% ...this happens, e.g., in @nonlinearSys with reach (from contDynamics)

% safe overview
save('unitTestsStatus.mat','overview');

% return to original working directory
cd(currentDirectory);

% execution successful
res = true;

end


% Auxiliary Functions -----------------------------------------------------

function funclist = functionList(funcs)
% funcs must be struct with fields "dir", "files"
% caution: recursive calls

funclist = {};

for i=1:length(funcs.files)
    if isstruct(funcs.files{i}) % directory / class
        temp = functionList(funcs.files{i});
        funclist = [funclist; temp];
    else % file
        funclist{end+1,1} = [funcs.dir filesep funcs.files{i}];
    end
end


end

function unittestname = findUnitTest(func,testlist)
% find unit tests for given function (cell-array), admissible syntax:
%   test_classname_functionname(_...)   if function part of a class
%   test_functionname(_...)             if function not part of a class

% default: no match
unittestname = '---';

% read out function name (after last filesep)
idxfilesep = strfind(func,filesep);
funcname = extractBetween(func,idxfilesep(end)+1,length(func));
funcname = funcname{1};

% check whether function part of a class or not
if contains(func,'@')
    % read out class name
    idxClassStart = strfind(func,'@') + 1;
    idxClassEnd = idxfilesep(find(idxfilesep > idxClassStart,1,'first')) - 1;
    classname = extractBetween(func,idxClassStart,idxClassEnd);
    classname = classname{1};
else
    classname = [];
end

% define required matching syntax
if ~isempty(classname)
    testname = ['test_' classname '_' funcname];
else
    testname = ['test_' funcname];
end

temp = {};
% go through list to search for all matching names
for i=1:length(testlist)
    % either names are identical (shortest version), cases like
    %       test_linearSys_reachInner
    % or there is additional text (separated by a '_') as in
    %       test_linearSys_reach_01
    if strcmp(testlist{i},testname) || contains(testlist{i},[testname '_'])
        temp = [temp testlist{i}];
    end
end
if ~isempty(temp)
    unittestname = temp;
end


end

%------------- END OF CODE --------------
