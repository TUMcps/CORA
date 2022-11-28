function runTestSuite_nn(varargin)
% runTestSuite_nn - runs the standard test suite by executing all functions
% starting with the prefix 'test_nn'
%
% Syntax:  
%    runTestSuite_nn(varargin)
%
% Inputs:
%    verbose - show workspace output or not (not required)
%    directory - change directory (not required)
%
% Outputs:
%    -

% Author:       Tobias Ladner
% Written:      03-March-2022
% Last update:  28-November-2022 (testnn)
% Last revision:---


%------------- BEGIN CODE --------------

directory = [CORAROOT filesep 'unitTests'];
verbose = 1;

if nargin >= 1
    directory = varargin{1};
end

if nargin >= 2
    verbose = varargin{2};
end

% run main program performing the tests
% all 'neuralNetwork' nn tests  without additional toolboxes
[failed1, numberOfTests1] = testSuiteCore('test_nn',verbose,directory);
% all 'neurNetContrSys' nn tests without additional toolboxes
[failed2, numberOfTests2] = testSuiteCore('test_neurNetContrSys',verbose,directory);
% all tests requiring additional toolboxes
[failed3, numberOfTests3] = testSuiteCore('testnn',verbose,directory);

failed = [failed1; failed2; failed3];
numberOfTests = numberOfTests1 + numberOfTests2 + numberOfTests3;

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed, 1)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------