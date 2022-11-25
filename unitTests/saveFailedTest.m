function saveFailedTest (varargin)
% saveFailedTest - Saves the relevant data of a failed test so the occured
% error can be reviewed. This function should be called inside a test
% function once the failing of the test is apparant.
%
% Note: Function is primarily targeted at tests with random parameters.
%
%
% Syntax:  
%    saveFailedTest(varargin)
%
% Inputs:
%
%   varargin - Any data relevant to recreate the failed test-setup
%
% Outputs:
%    (.mat file with all relevant data, path of failed test and datetime of failed test
%     saved to unitTests/failedTests)
%
% Example: 
%    -

% Author:       Maximilian Perschl
% Written:      02-August-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% get file name of test function
stack = dbstack(1);
testFileName = stack.name;

% get current date and time as string
timestr = datestr(datetime);
% change timestring to be compatible with filename-format
timestr = strrep(timestr,':','-');

% create .mat filename
filename = strcat(coraroot,'/unitTests/failedTests/',testFileName,'-',timestr,'.mat');

% save .mat file with all input arguments
save(filename,'varargin');


end