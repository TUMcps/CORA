function printTestOverview(varargin)
% printTestOverview - print information about the current state
%    of unit testing to the console
%
% Syntax:  
%    printTestOverview(varargin)
%
% Inputs:
%    format - (optional) what should be displayed:
%               'lastRun': results of last full test suite run
%               'classes': tests ordered by classes
%
% Outputs:
%    -
%
% Example: 
%    -

% Author:       Mark Wetzlinger
% Written:      22-January-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default format: last run
format = 'lastrun';

% parse input
allowedFormats = {'lastrun'; 'classes'};
if nargin == 1
    if ischar(varargin{1}) && any(ismember(allowedFormats,varargin{1}))
        format = varargin{1};
    end
end

% load data
if ~isfile([coraroot filesep 'unitTests' filesep 'unitTestsStatus.mat'])
    error("No data provided. Run '''runTestSuite''' to acquire data.");
else
    data = load('unitTestsStatus.mat');
    testdata = data.overview;
end


% output to console
if strcmp(format,'lastrun')
    
    % print header
    fprintf('-*---------------------------------*-\n');
    fprintf('--- Results of last full test run ---\n\n');

    % print machine and Matlab version
    fprintf(['  system:    ',testdata.lastRun.hostsystem,'\n']);
    fprintf(['  Matlab:    ',testdata.lastRun.matlabversion,'\n']);
    
    % print date and time
    fprintf(['  date/time: ',testdata.lastRun.date,'\n']);

    % print test results
    fprintf('  tests:\n');
    fprintf(['  .. total:  ',num2str(testdata.lastRun.nrTotalTests),'\n']);
    fprintf(['  .. failed: ',num2str(testdata.lastRun.nrFailedTests),'\n']);

    % print names (in case there are failed tests)
    if testdata.lastRun.nrFailedTests > 0
        fprintf(['  .. .. ',testdata.lastRun.nameFailedTests{1},'\n']);
        for i=2:testdata.lastRun.nrFailedTests
        fprintf(['        ',testdata.lastRun.nameFailedTests{i},'\n']);
        end
    end

    % print footer
    fprintf('-*---------------------------------*-\n');
    
elseif strcmp(format,'classes')
    % only goes through fullList
    
    % get filesep from system which produced the last run
    if contains(testdata.fullList.func{1},'/')
        delim = '/';
    elseif contains(testdata.fullList.func{1},'\')
        delim = '\';
    else
        error("Issue with last test run.");
    end        
    
    % truncate lists so that only functions within classes remain
    removeIdx = false(length(testdata.fullList.func),1);
    for i=1:length(testdata.fullList.func)
        if ~contains(testdata.fullList.func{i},'@')
            removeIdx(i) = true;
        end
    end
    funcs = testdata.fullList.func(~removeIdx);
    tests = testdata.fullList.test(~removeIdx);
    status = testdata.fullList.status(~removeIdx);
    
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

end

%------------- END OF CODE --------------