function runTestSuite_longDuration( varargin )
% runTestSuite_longDuration - runs the test suite with test that have a 
% long duration; all functions starting with the prefix 'testLongDuration_' 
% are executed
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

directory = [coraroot '/unitTests'];
verbose = 1;

if nargin >= 1
    directory = varargin{1};
end

if nargin >= 2
    verbose = varargin{2};
end

% run main program performing the tests
[failed, numberOfTests] = testSuiteCore('testLongDuration',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed, 2)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------