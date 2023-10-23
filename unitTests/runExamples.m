function runExamples(varargin)
% runExamples - runs all examples (files starting with 'example') in the
%    folder cora/examples
%
% Syntax:
%    runExamples(varargin)
%
% Inputs:
%    verbose - (optional) show workspace output or not
%    directory - (optional) change directory
%
% Outputs:
%    -

% Authors:       Dmitry Grebenyuk, Matthias Althoff
% Written:       20-September-2016
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get the original working directory
currentDirectory = pwd;

% default settings
[verbose,directory] = setDefaultValues(...
    {true,[CORAROOT filesep 'examples']},varargin);

% run main program performing the examples
[failed, numberOfTests] = testSuiteCore('example',verbose,directory);

% display number of failed examples and file names
disp('----------------------------------------------------------------------------');
disp(['run ' int2str(numberOfTests) ' examples, ' int2str(size(failed,1)) ' failed.']);
disp(strjoin(failed, ',\n'));

% return to original working directory
cd(currentDirectory);

% ------------------------------ END OF CODE ------------------------------
