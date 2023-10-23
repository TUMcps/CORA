function files = excludefiles(files,path,name)
% excludefiles - finds all files in a given path
%
% Syntax:
%    files = excludefiles(files,path,name)
%
% Inputs:
%    files - list of files
%    path - directory in which to exclude files
%    name - (optional) name of certain file to exclude
%
% Outputs:
%    files - updated list of files
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Tobias Ladner
% Written:       18-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

idx = contains({files.folder}', [CORAROOT filesep path]);
if nargin == 3
    idx = idx & contains({files.name}', name);
end
files = files(~idx);

% ------------------------------ END OF CODE ------------------------------
