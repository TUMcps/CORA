function runTestSuite(varargin)
% runTestSuite - runs the standard test suite by executing all functions
%    starting with the prefix 'test_'
%
% Syntax:  
%    runTestSuite(varargin)
%
% Inputs:
%    verbose - show workspace output or not (not required)
%    directory - change directory (not required)
%
% Outputs:
%    (text to console)
%
% Example: 
%    -

% Author:       Matthias Althoff, Mark Wetzlinger
% Written:      31-August-2016
% Last update:  22-January-2021 (MW, save results of full run)
% Last revision:---

%------------- BEGIN CODE --------------

% get the original working directory
currentDirectory = pwd;

% directory of unit tests
directory = [coraroot filesep 'unitTests'];
verbose = true;

if nargin >= 1
    verbose = varargin{1};
end
if nargin >= 2
    % correct separator
    temp = strrep(varargin{2},'\',filesep);
    temp = strrep(temp,'/',filesep);
    directory = temp;
    if ~isfolder(directory)
        error("Input given for directory is not a directory!");
    end
end

% save start date and time (only full runs)
if strcmp(directory,[coraroot filesep 'unitTests'])
    timeStart = char(datetime(now,'ConvertFrom','datenum'));
end

% run main program performing the tests
[failed, numberOfTests] = testSuiteCore('test',verbose,directory);

% save result (only full runs)
if strcmp(directory,[coraroot filesep 'unitTests'])
    timeEnd = char(datetime(now,'ConvertFrom','datenum'));
    % save end date and time, number (and names) of (failed) tests
    res_testrun.date = [timeStart ' - ' timeEnd(13:end)];
    res_testrun.nrTotalTests = numberOfTests;
    res_testrun.nrFailedTests = size(failed,1);
    res_testrun.nameFailedTests = failed;
    % save matlab version
    res_testrun.matlabversion = version;
    % save operating system
    res_testrun.hostsystem = 'Unknown';
    if ismac
        res_testrun.hostsystem = 'Mac';
    elseif isunix
        res_testrun.hostsystem = 'Unix';
    elseif ispc
        res_testrun.hostsystem = 'Windows';
    end
    % save result
    unitTestsOverview(res_testrun);
    % print result
    printTestOverview('lastrun');
    
else
    % display number of failed tests + file names
    disp('----------------------------------------------------------------------------');
    disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed,1)) ' failed.']);
    disp(strjoin(failed, ',\n'));
    
end

% return to original working directory
cd(currentDirectory);

%------------- END OF CODE --------------