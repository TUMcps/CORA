function runTestSuite_Mosek(varargin)
% runTestSuite_Mosek - runs the test suite with test that require the Mosek
% solver. The tests have been designed for the Mosek 2021 version; all 
% functions starting with the prefix 'testMosek_' are executed
%
% Syntax:  
%    runTestSuite_Mosek(varargin)
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
% Written:      07-July-2021
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
[failed, numberOfTests] = testSuiteCore('testMosek',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed, 2)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------