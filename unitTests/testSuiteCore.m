function [failed, numberOfTests] = testSuiteCore(varargin)
% testSuiteCore - runs functions starting with a certain prefix contained 
% in the directory dir or in a subdir
%
% Syntax:  
%    [failed, numberOfTests] = testSuiteCore(varargin)
%
% Inputs:
%    prefix - prefix of function names to be tested
%    verbose - show workspace output or not (not required)
%    directory - change directory (not required)
%
% Outputs:
%    failed - a cell array containing the name of the failed tests
%    numberOfTests - number of performed tests
%
% Example: 
%    -
%
% 
% Author:       Dmitry Grebenyuk, Matthias Althoff
% Written:      31-August-2016
% Last update:  ---
% Last revision:---


%------------- BEGIN CODE --------------

directory = [coraroot '/unitTests'];
verbose = 0;

if nargin >= 1
    prefix = varargin{1};
end
if nargin >= 2
    verbose = varargin{2};
end
if nargin >= 3
    directory = varargin{3};
end

numberOfTests = 0;
failed = cell(0);

% Find all relevant files in directory und run the tests
files = dir([directory, ['/',prefix,'_*.m']]);
for i=1:size(files,1)
    % Extract the function name
    [~, fname] = fileparts(files(i).name);
    % Supress output of tests by usage of evalc
    try
        fprintf(['run ',fname,': ']);
        [~,res] = evalc(fname);
        if res == 1
            fprintf('%s\n','passed');
        else
            fprintf('%s\n','failed');
        end
    catch
        res = 0;
        fprintf('%s\n','failed');
    end
    
    if res == 0
        failed = [failed {fname}];
    end
    numberOfTests = numberOfTests + 1;
end

% run files in subdirectories
files = dir(directory);
for i=1:size(files,1)
    % Exclude files and directories . and ..
    if files(i).name(1) == '.' || files(i).isdir == 0
        continue;
    end
        
    [subfailed, subnum] = testSuiteCore(prefix, verbose, [directory '/' files(i).name]);
    failed = [failed subfailed];
    numberOfTests = numberOfTests + subnum;
end

% if verbose ~= 0
%     disp('----------------------------------------------------------------------------');
%     disp(['run ' int2str(numberOfTests) ' tests, ' int2str(size(failed, 2)) ' failed.']);
%     disp(strjoin(failed, ',\n'));
% end

end

%------------- END OF CODE --------------