function testsByClasses
% testsByClasses - print on command window the number of tests for each
%    CORA class
%
% Syntax:
%    testsByClasses
%
% Inputs:
%    -
%
% Outputs:
%    -

% Authors:       Mark Wetzlinger
% Written:       09-April-2023 (moved from printTestOverview and unitTestsOverview)
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% switch to directory of unit tests
directory = [CORAROOT filesep 'unitTests'];
cd(directory);

% check if file available (as input argument or saved .mat file)
if ~isfile('unitTestsStatus.mat')
    throw(CORAerror('CORA:specialError',...
        "No data provided. Run '''runTestSuite''' to acquire data."));
else
    % load results from last full test run
    load('unitTestsStatus.mat','testResults');
    % results stored in map with key = 'short', value = struct with data
    testdata = testResults('short');
end

% list main content folders
foldernames = {'contDynamics';'contSet';'discrDynamics';'hybridDynamics';'matrixSet'};

% name of failed tests
nameFailedTests = {testdata.results.fname}';
nameFailedTests = nameFailedTests(~[testdata.results.ok]');

% write all funcs in one list and find matching unit tests
allfuncs = {};
alltests = {};
for i=1:length(foldernames)
    % list functions and unit tests
    tempfiles = listFolderContent([CORAROOT filesep foldernames{i}]);
    temptests = listFolderContent([CORAROOT filesep 'unitTests' filesep foldernames{i}]);
    % collapse lists into one list
    tempfiles = aux_functionList(tempfiles);
    allfuncs = [allfuncs; tempfiles];
    temptests = aux_functionList(temptests);
    alltests = [alltests; temptests];
end

corarootlength = length(CORAROOT);
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

fullList.func = allfuncs;
% find unit test for each function
for i=1:length(allfuncs)
    fullList.test{i,1} = aux_findUnitTest(allfuncs{i},alltests);
    % check if run was successful in last run
    if strcmp(fullList.test{i},'---')
        fullList.status(i,1) = -1; % no information
    else
        if any(ismember(nameFailedTests,fullList.test{i}))
            fullList.status(i,1) = 0; % at least one test failed
        else
            fullList.status(i,1) = 1; % all tests successful
        end
        % delete test(s) from alltests
        for j=1:length(fullList.test{i})
            alltests = alltests(~strcmp(alltests,fullList.test{i}{j}));
        end
        fullList.test{i} = strjoin(fullList.test{i},'\n');
    end
end

% unassigned test: either duplicates (e.g. reach_02) or name not found in func
unassignedList.test = alltests;
if ~isempty(unassignedList.test)
unassignedList.status = ...
    ~ismember(unassignedList.test,nameFailedTests);
end
% issue: if method from superclass called, function test in unassignedList
% ...this happens, e.g., in @nonlinearSys with reach (from contDynamics)

% get filesep from system which produced the last run
if contains(fullList.func{1},'/')
    delim = '/';
elseif contains(fullList.func{1},'\')
    delim = '\';
else
    throw(CORAerror('CORA:specialError',...
        'Issue with last test run.'));
end        

% truncate lists so that only functions within classes remain
removeIdx = false(length(fullList.func),1);
for i=1:length(fullList.func)
    if ~contains(fullList.func{i},'@')
        removeIdx(i) = true;
    end
end
funcs = fullList.func(~removeIdx);
tests = fullList.test(~removeIdx);
status = fullList.status(~removeIdx);

% read out superfolders of classes (since no class on root level)
firstfilesep = strfind(funcs{1},delim);
firstfilesep = firstfilesep(1);
superfolders = extractBetween(funcs{1},1,firstfilesep-1);
currsuperfolder = superfolders{1}; startIdx = 1;
c = 1;
for i=1:length(funcs)
    firstfilesep = strfind(funcs{i},delim);
    firstfilesep = firstfilesep(1);
    temp = funcs{i}(1:firstfilesep-1); 
    if ~strcmp(currsuperfolder,temp)
        endIdx(c,1) = i-1; c = c+1; startIdx(c,1) = i;
        superfolders{c,1} = temp;
        currsuperfolder = superfolders{c};
    end
end
endIdx(c,1) = length(funcs);

currclass = '';
% read classes / nr of functions per class inside superfolders
for i=1:length(superfolders)
    c = 0; nrfunc = 0;
    for j=startIdx(i):endIdx(i)
        idxatsign = strfind(funcs{j},'@');
        idxatsign = idxatsign(1);
        idxfileseps = strfind(funcs{j},delim);
        nextfilesep = idxfileseps(find(idxfileseps > idxatsign,1,'first'));
        temp = funcs{j}(idxatsign+1:nextfilesep-1);
        if ~strcmp(currclass,temp)
            if c > 0; nrfuncperclass{i,1}(c,1) = nrfunc; end
            c = c+1;
            classes{i,1}{c} = temp;
            currclass = temp;
            nrfunc = 0;
        end
        nrfunc = nrfunc + 1;
    end
    nrfuncperclass{i}(c,1) = nrfunc;
end

% read number of tests / failed tests
for i=1:length(superfolders)
    idx = startIdx(i);
    for c=1:length(nrfuncperclass{i})
        nrtestsperclass{i,1}(c,1) = nnz(status(idx:idx+nrfuncperclass{i}(c)-1) ~= -1);
        nrfailedtestsperclass{i,1}(c,1) = nnz(status(idx:idx+nrfuncperclass{i}(c)-1) == 0);
        idx = idx + nrfuncperclass{i}(c);
    end
    nrtestspersuperfolder{i,1} = num2str(sum(nrtestsperclass{i}));
    nrfailedtestspersuperfolder{i,1} = num2str(sum(nrfailedtestsperclass{i}));
end

% find out max lengths of char arrays for table
maxSuperfolder = 0;
maxClass = 0;
for i=1:length(superfolders)
    if length(superfolders{i}) > maxSuperfolder
        maxSuperfolder = length(superfolders{i});
    end
    for j=1:length(classes{i})
        if length(classes{i}{j}) > maxClass
            maxClass = length(classes{i}{j});
        end
    end
end


% print header
fprintf('-*---------------------------------------------------*-\n');
fprintf('---              Unit tests by classes              ---\n\n');

% count char arrays for vertical alignment
indent = max(maxSuperfolder,maxClass+3);
nrfunc = 'nr.func  ';
nrtests = 'nr.tests  ';
nrfailed = 'nr.failed';
nrfuncblock = length(nrfunc);
nrtestsblock = length(nrtests);

% first line
fprintf([repmat(' ',1,indent) ' ' nrfunc nrtests nrfailed '\n']);

% loop over superfolders and classes
for i=1:length(superfolders)
    % print name of superfolder
    fprintf([superfolders{i} repmat(' ',1,indent-length(superfolders{i})+1)]);
    nrFuncFolder = num2str(endIdx(i)-startIdx(i)+1);
    % print number of functions
    fprintf([nrFuncFolder repmat(' ',1,nrfuncblock-length(nrFuncFolder))]);
    % print number of tests
    fprintf([nrtestspersuperfolder{i} ...
        repmat(' ',1,nrtestsblock-length(nrtestspersuperfolder{i}))]);
    % print number of failed tests
    fprintf([nrfailedtestspersuperfolder{i} '\n']);
    
    for j=1:length(classes{i})
        % print name of class
        fprintf(['.. ' classes{i}{j} repmat(' ',1,indent-length(classes{i}{j})-2)]);
        % print number of functions
        fprintf([num2str(nrfuncperclass{i}(j)) ...
            repmat(' ',1,nrfuncblock-length(num2str(nrfuncperclass{i}(j))))]);
        % print number of tests
        fprintf([num2str(nrtestsperclass{i}(j)) ...
            repmat(' ',1,nrtestsblock-length(num2str(nrtestsperclass{i}(j))))]);
        % print number of failed tests
        fprintf([num2str(nrfailedtestsperclass{i}(j)) '\n']);
    end
    
end

% print footer
fprintf('-*---------------------------------------------------*-\n');

end


% Auxiliary functions -----------------------------------------------------

function funclist = aux_functionList(funcs)
% funcs must be struct with fields "dir", "files"
% caution: recursive calls

funclist = {};
for i=1:length(funcs.files)
    if isstruct(funcs.files{i}) % directory / class
        temp = aux_functionList(funcs.files{i});
        funclist = [funclist; temp];
    else % file
        funclist{end+1,1} = [funcs.dir filesep funcs.files{i}];
    end
end

end

function unittestname = aux_findUnitTest(func,testlist)
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

% ------------------------------ END OF CODE ------------------------------
