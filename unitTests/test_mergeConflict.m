function res = test_mergeConflict
% test_mergeConflict - check if any merge conflicts were left unresolved
%
% Syntax:
%    res = test_mergeConflict
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_docstring

% Authors:       Mark Wetzlinger
% Written:       28-April-2023 
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
    
% get all files
files = findfiles(CORAROOT);
% exclude this file (written a bit cumbersome to deal with future moves)
thisFolder = strrep(mfilename('fullpath'),[CORAROOT filesep],'');
lastFilesep = strfind(thisFolder,filesep);
lastFilesep = lastFilesep(end);
thisFolder = thisFolder(1:lastFilesep-1);
files = excludefiles(files,thisFolder,[mfilename '.m']);

% list unresolved files
unresolvedFiles = {};

% loop over all files
for i=1:length(files)
    file = files(i);
    file_path = [file.folder, filesep, file.name];
    filetext = fileread(file_path);

    % possibly rewrite the check below in case these strings turn up
    % somewhere else unrelated to merge conflict resolutions
    if contains(filetext,">>>>>>>") ...
            || contains(filetext,"<<<<<<< HEAD:")
        % add to list
        unresolvedFiles{end+1} = file.name;

        % display information
        disp("Suspected unresolved merge conflict in " ...
            + file.name);
    end
end

% test is ok if no remaining unresolved merge conflicts
res = isempty(unresolvedFiles);

% ------------------------------ END OF CODE ------------------------------
