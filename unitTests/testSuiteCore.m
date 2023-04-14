function [failed, numberOfTests] = testSuiteCore(prefix,varargin)
% testSuiteCore - runs functions starting with a certain prefix contained 
%    in the directory and recursively searches all subfolders for same prefix
%
% Syntax:  
%    [failed, numberOfTests] = testSuiteCore(prefix)
%    [failed, numberOfTests] = testSuiteCore(prefix,verbose)
%    [failed, numberOfTests] = testSuiteCore(prefix,verbose,directory)
%
% Inputs:
%    prefix - prefix of function names to be tested
%    verbose - (optional) show workspace output or not
%    directory - (optional) change directory
%
% Outputs:
%    failed - a cell array containing the names of the failed tests
%    numberOfTests - number of run tests
%
% Example:
%    -

% Author:       Dmitry Grebenyuk, Matthias Althoff
% Written:      31-August-2016
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default values
[verbose,directory] = setDefaultValues(...
    {false,[CORAROOT filesep 'unitTests']},varargin);

numberOfTests = 0;
failed = cell(0);

% Find all relevant files in directory und run the tests
files = dir([directory filesep prefix '_*.m']);
for i=1:size(files,1)
    % Extract the function name
    [~, fname] = fileparts(files(i).name);
    % Supress output of tests by usage of evalc
    try
        fprintf(['run ' fname ': ']);
        [~,res] = evalc(fname);
        if res
            fprintf('%s\n','passed');
        else
            fprintf('%s\n','failed');
        end
    catch
        res = false;
        fprintf('%s\n','failed');
    end
    
    % save file names of failed tests
    if ~res
        failed = [failed; {fname}];
    end
    numberOfTests = numberOfTests + 1;
end

% run files in subdirectories
files = dir(directory);
for i=1:size(files,1)
    % Exclude directories "." and "..", as well as files (checked above)
    if files(i).name(1) == '.' || files(i).isdir == 0
        continue;
    end
    
    % recursive call of subfolder (and its subfolders, ... etc.)
    [subfailed, subnum] = testSuiteCore(prefix,verbose,[directory filesep files(i).name]);
    failed = [failed; subfailed];
    numberOfTests = numberOfTests + subnum;
end

%------------- END OF CODE --------------