function tree = listFolderContent(rootDir)
% listFolderContent - recursive function to list all functions in rootDir
%
% Syntax:
%    tree = listFolderContent(rootDir)
%
% Inputs:
%    rootDir - directory whose content is to be listed
%
% Outputs:
%    tree - tree structure with content of rootDir
%
% Example: 
%    -

% Authors:       Mark Wetzlinger
% Written:       20-January-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% files in rootDir
files = dir(rootDir);
% go through directories / files
fnames.dir = rootDir;
fnames.files = {};
for i=1:length(files)
    % exclude directories "." and ".." and private
    if strcmp(files(i).name(1),'.') || strcmp(files(i).name,'private')
        continue;
    end
    % regular directories (not private): call function again
    if files(i).isdir
        fnames.files{end+1,1} = ...
            listFolderContent([files(i).folder filesep files(i).name]);
        continue;
    end
    % extract the function name (only if .m-file)
    [~, fname,ext] = fileparts(files(i).name);
    if strcmp(ext,'.m')
        fnames.files{end+1,1} = fname;
    end
end

tree = fnames;

% ------------------------------ END OF CODE ------------------------------
