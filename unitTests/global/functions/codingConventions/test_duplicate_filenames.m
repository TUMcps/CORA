function res = test_duplicate_filenames()
% test_duplicate_filenames - tests for duplicate file names
%    these should be avoided for Matlab to call the correct function
%
% Syntax:
%    res = test_duplicate_filenames()
%
% Inputs:
%    -
%
% Outputs:
%    res - whether all files are valid
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       26-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% find all files
files = findfiles(CORAROOT,true);

% exclude:
% - thirdparty files
files = excludefiles(files, ['global' filesep 'thirdparty']);

% find unique filenames
filenames = unique({files.name});

% iterate through all file names
for i=1:numel(filenames)
    filename = filenames{i};

    % find files with that filename
    idx = ismember({files.name},filename);
    files_i = files(idx);

    % check if only one file was found
    if isscalar(files_i)
        continue
    end

    % check if all but one are part of a class
    numNonClassFiles = nnz(~contains({files_i.folder},'@'));
    if numNonClassFiles>1
        % found issue
        fprintf('Issue for ''%s'':\n',filename);
        for j=1:numel(files_i)
            fprintf('- <a href="matlab:open(''%s/%s'')">%s/%s</a>\n',files_i(j).folder,files_i(j).name,files_i(j).folder,files_i(j).name)
        end
        fprintf('\n')
        res = false;
    end
end

end

% ------------------------------ END OF CODE ------------------------------
