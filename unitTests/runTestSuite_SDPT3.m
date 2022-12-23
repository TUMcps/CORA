function runTestSuite_SDPT3(varargin)
% runTestSuite_SDPT3 - runs the test suite with test that require the SDPT3
%    solver. The tests have been designed for the SDPT3 V.4.0; all 
%    functions starting with the prefix 'testSDPT3_' are executed
%
% Syntax:  
%    runTestSuite_SDPT3(varargin)
%
% Inputs:
%    verbose - show workspace output or not (not required)
%    directory - change directory (not required)
%
% Outputs:
%    -

% Author:       Matthias Althoff
% Written:      07-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% too many input arguments
if nargin > 2
    throw(CORAerror('CORA:tooManyInputArgs',2));
end

% set default values
[directory,verbose] = setDefaultValues({[CORAROOT filesep 'unitTests'],true},varargin);


% run main program performing the tests
[failed, numberOfTests] = testSuiteCore('testSDPT3',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed,1)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------