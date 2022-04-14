function path = pathFailedTests(file)
% pathFailedTest - returns a filepath to the directory that stores all
%                  failed unit tests
%
% Syntax:  
%    path = pathFailedTests(file)
%
% Inputs:
%    file - filename of the unit-test that failed
%
% Outputs:
%    path - path to the corresponding file stroing data for the unit test
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: coraroot

% Author:       Niklas Kochdumper
% Written:      02-February-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------
    
    % create file name
    file_name = strcat(file,'_',datestr(now,'mm-dd-yyyy_HH-MM'));
    
    % create full path to the file
    path = fullfile(coraroot(), 'unitTests', 'failedTests',file_name);
    
    % create directory if it does not exist yet
    dirPath = fullfile(coraroot(),'unitTests','failedTests');
    
    if ~isdir(dirPath)
       mkdir(dirPath); 
    end
end

%------------- END OF CODE --------------