function res = test_docstring()
% test_docstring - tests if all CORA files have a valid docstring
%
% Syntax:
%    res = test_docstring()
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

% Author:        Tobias Ladner
% Written:       18-November-2022
% Last update:   ---
% Last revision: ---

%------------- BEGIN CODE --------------

files = [
    aux_addFiles('contDynamics');
    aux_addFiles('contSet');
    aux_addFiles('converter');
    aux_addFiles('contDynamics');
    % aux_addFiles('discrDynamics');
    aux_addFiles('examples');
    aux_addFiles('global');
    aux_addFiles('hybridDynamics');
    aux_addFiles('matrixSet');
    aux_addFiles('unitTests');
];

% exclude file paths
% converter
files = aux_exclFiles(files, ['converter' filesep 'powerSystem2cora' filesep 'cases'], 'IEEE14Parameters.m');
files = aux_exclFiles(files, ['converter' filesep 'powerSystem2cora' filesep 'cases'], 'IEEE30Parameters.m');
% examples
files = aux_exclFiles(files, ['examples' filesep 'manual']);
% global
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'combinator']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'cprnd']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'eq_sphere_partitions']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'm2tex']);
files = aux_exclFiles(files, ['global' filesep 'functions' filesep 'tprod']);
files = aux_exclFiles(files, ['global' filesep 'functions'], 'Direct.m');
% unitTests
files = aux_exclFiles(files, ['unitTests' filesep 'converter' filesep 'powerSystem2cora' filesep 'models']);

% check in reverse order
% files = flipud(files);

res = true;
issue_count = 0;
for i=1:length(files)
    file = files(i);
    file_path = [file.folder, filesep, file.name];
    filetext = fileread(file_path);

    issues = [];

    % check docstring

    if ~contains(filetext, "BEGIN CODE")
        issues{end+1} = "'BEGIN CODE' missing";
    end

    if ~contains(filetext, "END OF CODE")
        issues{end+1} = "'END OF CODE' missing";
    end

    if count(filetext, file.name(1:end-2)) < 3
        issues{end+1} = "Filename should appear at least 3 times (function/class name, begin of docstring, syntax)";
    end

    % if filetext(1:8) == 'classdef'
    %     if ~contains(filetext, "% Description")
    %        issues{end+1} = "Description missing";
    %     end
    % end

    if ~contains(filetext, "% Syntax")
        issues{end+1} = "Syntax missing.";
    end

    if ~contains(filetext, "% Input")
        issues{end+1} = "Input missing.";
    end

    if ~contains(filetext, "% Output")
        issues{end+1} = "Output missing.";
    end

    if ~contains(filetext, "% Author")
        issues{end+1} = "Author missing.";
    end

    authorPos = strfind(filetext, "Author");
    if ~isempty(authorPos)
        authorPos = authorPos(1);
        if filetext(authorPos-4) ~= filetext(authorPos-6)
            issues{end+1} = "Author should not be included in help text.";
        end
    end

    if ~contains(filetext, "% Written")
        issues{end+1} = "Written missing.";
    end

    if ~contains(filetext, "% Last update")
        issues{end+1} = "Last update missing.";
    end

    if ~contains(filetext, "% Last revision")
        issues{end+1} = "Last revision missing.";
    end

    % display information

    all_good = isempty(issues);
    if ~all_good
        disp(['- ' file_path])
        disp("  Issues:")
        fprintf("  - %s\n", string(issues'))
        issue_count = issue_count + 1;
        res = false;
        % might be helpful to directly open wrong file in the debugger
        % open(file_path)
        disp(" ")
    end
end

if issue_count > 0
    fprintf("Found errors in docstring in %i files.\n", issue_count);
end

end

% Auxiliary Functions -----------------------------------------------------

function files = aux_addFiles(path, includeSubfolders)
    if nargin < 2
        includeSubfolders = true;
    end

    subpath = '';
    if includeSubfolders
        subpath = ['**' filesep];
    end

    files = dir([CORAROOT filesep path filesep subpath '*.m']);
end

function files = aux_exclFiles(files, path, name)
    idx = contains({files.folder}', [CORAROOT filesep path]);
    if nargin == 3
        idx = idx & contains({files.name}', name);
    end
    files = files(~idx);
end

%------------- END OF CODE --------------