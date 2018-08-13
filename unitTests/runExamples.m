function runExamples(varargin)
% runExamples - runs all examples starting with a certain prefix in the
% folder 'examples'.
%
% Syntax:  
%    runExamples(varargin)
%
% Inputs:
%    prefix - prefix of function names to be tested
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
% Author:       Dmitry Grebenyuk, Matthias Althoff
% Written:      20-September-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

directory = [coraroot '/examples'];
prefix = 'example';
verbose = 1;

if nargin >= 1
    prefix = varargin{1};
end
if nargin >= 2
    verbose = varargin{2};
end
if nargin >= 3
    directory = varargin{3};
end

% run main program performing the tests
[failed, numberOfTests] = testSuiteCore(prefix,verbose,directory);

disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' examples, ' int2str(size(failed, 2)) ' failed.']);
disp(strjoin(failed, ',\n'));

%return to original working directory
cd(directory);

%------------- END OF CODE --------------