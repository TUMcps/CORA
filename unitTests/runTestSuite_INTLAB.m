function runTestSuite_INTLAB( varargin )
% runTestSuite_INTLAB - runs the standard test suite by executing all functions
% starting with the prefix 'testINTLAB_'
%
% Syntax:  
%    runTestSuite_INTLAB( varargin )
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
[failed, numberOfTests] = testSuiteCore('testINTLAB',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed, 2)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------