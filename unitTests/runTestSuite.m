function runTestSuite( varargin )
% runTestSuite - runs the standard test suite by executing all functions
% starting with the prefix 'test_'
%
% Syntax:  
%    [failed, numberOfTests] = testSuiteCore(varargin)
%
% Inputs:
%    verbose - show workspace output or not (not required)
%    directory - change directory (not required)
%
% Outputs:
%    -
%
% Example: 
%    -
%
% 
% Author:       Matthias Althoff
% Written:      31-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

% get the original working directory
currentDirectory = pwd;

directory = [coraroot '/unitTests'];
verbose = 1;

if nargin >= 1
    directory = varargin{1};
end

if nargin >= 2
    verbose = varargin{2};
end

% run main program performing the tests
[failed, numberOfTests] = testSuiteCore('test',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed, 2)) ' failed.']);
disp(strjoin(failed, ',\n'));

%return to original working directory
cd(currentDirectory);

%------------- END OF CODE --------------