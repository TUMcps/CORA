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
% Last update:   21-April-2023 (unix bugfix: author follows empty line)
%                19-July-2023 (ignore script files)
% Last revision: ---

%------------- BEGIN CODE --------------

files = [
    findfiles([CORAROOT filesep 'contDynamics']);
    findfiles([CORAROOT filesep 'contSet']);
    findfiles([CORAROOT filesep 'converter']);
    findfiles([CORAROOT filesep 'contDynamics']);
    % findfiles([CORAROOT filesep 'discrDynamics']);
    findfiles([CORAROOT filesep 'examples']);
    findfiles([CORAROOT filesep 'global']);
    findfiles([CORAROOT filesep 'hybridDynamics']);
    findfiles([CORAROOT filesep 'matrixSet']);
    findfiles([CORAROOT filesep 'unitTests']);
];

% exclude file paths
% converter
files = excludefiles(files, ['converter' filesep 'powerSystem2cora' filesep 'cases'], 'IEEE14Parameters.m');
files = excludefiles(files, ['converter' filesep 'powerSystem2cora' filesep 'cases'], 'IEEE30Parameters.m');
% examples
files = excludefiles(files, ['examples' filesep 'manual']);
% thirdparty
files = excludefiles(files, ['global' filesep 'thirdparty']);
% unitTests
files = excludefiles(files, ['unitTests' filesep 'converter' filesep 'powerSystem2cora' filesep 'models']);

% check in reverse order
% files = flipud(files);

res = true;
issue_count = 0;
for i=1:length(files)
    file = files(i);
    file_path = [file.folder, filesep, file.name];
    filetext = fileread(file_path);

    if ~startsWith(filetext,'function') && ~startsWith(filetext,'classdef')
        % ignore script files
        continue
    end

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
    else
        % author text should follow empty line
        authorPos = strfind(filetext, "% Author");
        lines = split(filetext(1:authorPos(1)), compose("\r\n"));
        if length(lines) == 1
            lines = split(filetext(1:authorPos(1)), newline);
        end
        if ~isempty(lines{end-1})
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

    if startsWith(file.folder, [CORAROOT filesep 'unitTests' filesep])
        if contains(filetext, 'close all')
            issues{end+1} = "Avoid 'close all' in unit tests.";
        end
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

%------------- END OF CODE --------------
