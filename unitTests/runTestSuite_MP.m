function runTestSuite_MP(varargin)
% runTestSuite_MP - runs the test suite with test that require the 
%    multiple precision toolbox; all functions starting with the prefix 
%    'testMP_' are executed
%
% Syntax:  
%    runTestSuite_MP(varargin)
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
[failed, numberOfTests] = testSuiteCore('testMP',verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed,1)) ' failed.']);
disp(strjoin(failed, ',\n'));

%------------- END OF CODE --------------