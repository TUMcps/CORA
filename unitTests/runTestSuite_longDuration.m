function runTestSuite_longDuration(varargin)
% runTestSuite_longDuration - runs the test suite with test that have a 
% long duration; all functions starting with the prefix 'testLongDuration_' 
% are executed
%
% Syntax:  
%    runTestSuite_longDuration(varargin)
%
% Inputs:
%    verbose - show workspace output or not (not required)
%    directory - change directory (not required)
%
% Outputs:
%    -

% Author:       Matthias Althoff
% Written:      31-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% too many input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% set default values
[directory,verbose] = setDefaultValues({[CORAROOT filesep 'unitTests'],true},varargin{:});

% run main program performing the tests
[failed, numberOfTests] = testSuiteCore('testLongDuration',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed,1)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------