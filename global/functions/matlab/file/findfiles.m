function files = findfiles(path,varargin)
% findfiles - finds all files in a given path
%
% Syntax:
%    files = findfiles(path,includeSubfolders,prefix)
%
% Inputs:
%    path - directory in which to search for files
%    includeSubfolders - (optional) true/false whether to include subfolders
%    prefix - (optional) prefix for file name
%
% Outputs:
%    files - list of found files
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

% default values
[includeSubfolders,prefix] = setDefaultValues({true,''},varargin);

% default: no subpath
subpath = '';
if includeSubfolders
    % '**' searches all subfolders
    subpath = ['**' filesep];
end

% list files
files = dir([path filesep subpath prefix '*.m']);

% ------------------------------ END OF CODE ------------------------------
