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
% See also: fix_docstring

% Authors:       Tobias Ladner
% Written:       18-November-2022
% Last update:   21-April-2023 (unix bugfix: author follows empty line)
%                19-July-2023 (ignore script files)
%                19-January-2024 (stricter syntax check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
% files = excludefiles(files, ['examples' filesep 'manual']);
% thirdparty
files = excludefiles(files, ['global' filesep 'thirdparty']);
% unitTests
files = excludefiles(files, ['unitTests' filesep 'converter' filesep 'powerSystem2cora' filesep 'models']);

% check in reverse order
% files = flipud(files);

res = true;
issue_count = 0;
for i=1:length(files)
    
    % read file
    file = files(i);
    file_path = [file.folder, filesep, file.name];
    filename = file.name(1:end-2); 

    % read content of file
    filetext = fileread(file_path);

    % collect issues
    issues = {};

    % split lines
    lines = splitlines(filetext);

    % hack to use break within code checks
    everythingTested = true;
    for h=1 

        % initial checks --------------------------------------------------

        if isempty(lines)
            issues{end+1} = 'Empty file.';
            break
        end

        if ~startsWith(lines{1},'function') && ~startsWith(filetext,'classdef')
            % ignore script files
            continue
        end

        % init line counter
        lcnt = 1;

        try
            cnt = 0;
            while  lcnt <= length(lines) && ~startsWith(lines{lcnt},'%') && ...
                (contains(lines{lcnt},filename) || contains(lines{lcnt},'...') || contains(lines{lcnt-1},'...'))
                if contains(lines{lcnt},filename)
                    cnt = cnt + 1;
                end
                lcnt = lcnt + 1;
            end
        catch ME
            cnt = 0;
        end

        if cnt == 0
            issues{end+1} = "Filename and class/function are called differently.";
            break
        end

        if length(lines) < lcnt || ~startsWith(lines{lcnt},['% ' filename])
            issues{end+1} = "Docstring should start with '% <filename>'.";
        end
        lcnt = lcnt + 1;

        % test syntax block in docstring ----------------------------------

        % iterate to syntax block
        while lcnt <= length(lines) && ~strcmp(lines{lcnt}, '% Syntax:')
            lcnt = lcnt + 1;
        end

        if lcnt > length(lines)
            issues{end+1} = 'Syntax block is missing.';
            break
        end

        if ~isEmptyCommentLine(lines{lcnt-1})
            issues{end+1} = 'Empty line before syntax block is missing.';
        end
        lcnt = lcnt + 1;

        % should have at least one exemplary call
        cnt = 0;
        while lcnt <= length(lines) && startsWith(lines{lcnt},'%') && ~isEmptyCommentLine(lines{lcnt})
            if ~isempty(regexp(lines{lcnt}, sprintf('[\\.= ]%s(\\(.*|)$',filename), 'once','emptymatch'))
                cnt = cnt + 1;
            end
            lcnt = lcnt + 1;
        end
        if cnt == 0
            issues{end+1} = 'At least one exemplary call should be stated.';
        end

        if ~isEmptyCommentLine(lines{lcnt})
            issues{end+1} = 'Empty line after syntax block is missing.';
        end
        lcnt = lcnt + 1;

        % test inputs block in docstring ----------------------------------

        % iterate to output block
        while lcnt <= length(lines) && ~strcmp(lines{lcnt}, '% Inputs:')
            lcnt = lcnt + 1;
        end

        if lcnt > length(lines) || ~strcmp(lines{lcnt}, '% Inputs:')
            issues{end+1} = 'Inputs block is missing.';        
            break
        end        

        if ~isEmptyCommentLine(lines{lcnt-1})
            issues{end+1} = 'Empty line before inputs block is missing.';
        end
        lcnt = lcnt+1;

        if isEmptyCommentLine(lines{lcnt})
            issues{end+1} = 'Missing inputs description.';
        end
        lcnt = lcnt+1;

        
        % test outputs block in docstring ---------------------------------

        % iterate to output block
        while lcnt <= length(lines) && ~strcmp(lines{lcnt}, '% Outputs:')
            lcnt = lcnt + 1;
        end

        if lcnt > length(lines)
            issues{end+1} = 'Outputs block is missing.';        
            break
        end

        if ~isEmptyCommentLine(lines{lcnt-1})
            issues{end+1} = 'Empty line before outputs block is missing.';
        end
        lcnt = lcnt+1;

        if isEmptyCommentLine(lines{lcnt})
            issues{end+1} = 'Missing outputs description.';
        end
        lcnt = lcnt+1;

        % test authors block in docstring ---------------------------------

        % iterate to authors block
        while lcnt <= length(lines) && ~startsWith(lines{lcnt}, '% Author')
            lcnt = lcnt + 1;
        end

        if lcnt > length(lines)
            issues{end+1} = 'Authors block is missing.';        
            break
        end

        if ~startsWith(lines{lcnt}, '% Authors:       ')
            issues{end+1} = 'Authors block is written inconsistently.';        
        end

        % check if previous line is empty
        if ~isempty(lines{lcnt-1})
            issues{end+1} = "Authors block should not be included in help text. Add an empty line prior to the author block.";
        end
        lcnt = lcnt + 1;

        % Written
        if ~startsWith(lines{lcnt},'% Written:')
            issues{end+1} = "'Written' line is missing.";
            break
        end
        if ~startsWith(lines{lcnt},'% Written:       ')
            issues{end+1} = "'Written' line is written inconsistently.";
            break
        end
        if ~validateAuthorDateLine(lines{lcnt})
            issues{end+1} = "'Written' date/comment is written inconsistently.";
        end
        lcnt = lcnt + 1;

        % Last update
        if ~startsWith(lines{lcnt},'% Last update:')
            issues{end+1} = "'Last update' line is missing.";
            break
        end
        if ~startsWith(lines{lcnt},'% Last update:   ')
            issues{end+1} = "'Last update' line is written inconsistently.";
        end

        while lcnt <= length(lines) && ~startsWith(lines{lcnt},'% Last revision:')
            if ~validateAuthorDateLine(lines{lcnt})
                issues{end+1} = "'Last update' date/comment are written inconsistently.";
            end
            lcnt = lcnt + 1;
        end

        % Last revision
        if lcnt > length(lines)
            issues{end+1} = "'Last revision' line is missing.";
            break
        end
        if ~startsWith(lines{lcnt},'% Last revision: ')
            issues{end+1} = "'Last revision' line is written inconsistently.";
        end

        while lcnt <= length(lines) && ~isEmptyLine(lines{lcnt})
            if ~validateAuthorDateLine(lines{lcnt})
                issues{end+1} = "'Last revision' date/comment are written inconsistently.";
            end
            lcnt = lcnt + 1;
        end
        lcnt = lcnt + 1;

        % test BEGIN CODE -------------------------------------------------

        if lcnt > length(lines)
            issues{end+1} = 'BEGIN CODE is missing.';        
            break
        end

        if ~strcmp(lines{lcnt}, '% ------------------------------ BEGIN CODE -------------------------------') ...
                || ~isEmptyLine(lines{lcnt-1}) || ~isEmptyLine(lines{lcnt+1})
            issues{end+1} = 'BEGIN CODE is written inconsistently.';
        end

        % test END CODE -------------------------------------------------

        % iterate to END CODE
        while lcnt <= length(lines) && ~contains(lines{lcnt}, ['END OF ' 'CODE'])
            lcnt = lcnt + 1;
        end

        if lcnt > length(lines)
            issues{end+1} = ['END OF ' 'CODE is missing.'];
            break
        end

        if ~strcmp(lines{lcnt}, ['% ------------------------------ END OF ' 'CODE ------------------------------']) ...
                || ~isEmptyLine(lines{lcnt-1}) ...
                || numel(lines) == lcnt || ~isEmptyLine(lines{lcnt+1})
            issues{end+1} = 'END CODE is written inconsistently.';
        end
        lcnt = lcnt + 1;

        if lcnt ~= length(lines) || ~isEmptyLine(lines{lcnt})
            issues{end+1} = ['END OF ' 'CODE is not at the end of file'];
        end

        % additional checks -----------------------------------------------
     
        % auxiliary functions
        lcnt = 2;
        while lcnt <= length(lines) 
            if startsWith(lines{lcnt}, 'function ')
                % check if it is a nested function
                if ~contains(lines{lcnt},'nest_')
                    % first auxiliary function found
                    break
                end
            end
            lcnt = lcnt + 1;
        end

        if lcnt <= length(lines)
            % first auxiliary function found, check previous lines
            
            while lcnt > 0 && ~strcmp(lines{lcnt}, '% Auxiliary functions -----------------------------------------------------')
                lcnt = lcnt - 1;
            end

            if lcnt == 0
                issues{end+1} = 'Start of auxiliary functions is missing.';
                break
            end

            if ~isEmptyLine(lines{lcnt-1}) || ~isEmptyLine(lines{lcnt-2}) || ~isEmptyLine(lines{lcnt+1}) || ...
                ~strcmp(lines{lcnt},'% Auxiliary functions -----------------------------------------------------')
                issues{end+1} = "Start of auxiliary functions is written inconsistently.";
            end        
            
            % auxiliary functions should start with 'aux_'
            while lcnt <= length(lines)
                if startsWith(lines{lcnt}, 'function ')

                    while ~contains(lines{lcnt},'aux_') && contains(lines{lcnt},'...')
                        lcnt = lcnt + 1;
                    end

                    if ~contains(lines{lcnt},'aux_')
                        issues{end+1} = "Auxiliary functions should start with 'aux_'.";
                        break
                    end
                end
                lcnt = lcnt + 1;
            end    
        end

        % avoid 'close all' in unit tests
        if startsWith(file.folder, [CORAROOT filesep 'unitTests' filesep])
            if contains(filetext, 'close all')
                issues{end+1} = "Avoid 'close all' in unit tests.";
            end
        end

        % avoid more than 2 empty lines
        cnt = 0; MAX_EMPTY = 2;
        for lcnt = 1:length(lines)
            if isEmptyCommentLine(lines{lcnt}) || isEmptyLine(lines{lcnt})
                cnt = cnt + 1;
            else
                if cnt > MAX_EMPTY
                    issues{end+1} = "Avoid more than two empty lines.";
                    break
                end

                % reset counter
                cnt = 0;
            end
        end

        % check see also for function_ files

        if endsWith(file.name, '_.m')
            lcnt = 1;
            while lcnt <= length(lines) && ~startsWith(lines{lcnt}, '% See also')
                lcnt = lcnt + 1;
            end
    
            if lcnt > length(lines)
                issues{end+1} = "Functions ending with '_' should contain a 'See also' line referencing the respective function in contSet.";
            else
                if ~startsWith(lines{lcnt}, sprintf('%% See also: contSet/%s', file.name(1:(length(file.name)-3))))
                    issues{end+1} = "Functions ending with '_' should reference the respective function in contSet in their 'See also' line.";
                end
            end
        end

        % everything tested -----------------------------------------------

        everythingTested = true;

    end

    % check if everything was tested
    if ~everythingTested
        issues{end+1} = '... and potentially some others.';
    end

    % spelling

    if contains(lower(filetext),['un' 'kown'])
        issues{end+1} = "Misspelled 'unknown'.";
    end

    % display information
    all_good = isempty(issues);
    if ~all_good
        disp(['- ' file_path])
        disp("  Issues:")
        fprintf("  - %s\n", string(issues'))
        disp(" ")

        % update result
        issue_count = issue_count + 1;
        res = false;

        % helpful to directly open wrong file in the debugger
        % open(file_path)
        % fix_docstring(file_path)
        % keyboard % acts as breakpoint   
    end
end

if issue_count > 0
    fprintf("Found errors in docstring in %i files.\n", issue_count);
    fprintf("Try to execute: fix_docstring(<file_path>)\n")
end

end

% ------------------------------ END OF CODE ------------------------------
