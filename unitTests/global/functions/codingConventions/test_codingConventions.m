function res = test_codingConventions()
% test_codingConventions - tests if all CORA files follow the coding conventions:
%    https://gitlab.lrz.de/cps/cora/-/wikis/home/Coding-conventions
%
% Syntax:
%    res = test_codingConventions()
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
%                15-February-2024 (FL, also check specification directory)
%                10-May-2024 (improved usability in command window)
%                16-July-2024 (TL, checks for CORAwarning and CORAlinprog)
%                05-February-2025 (TL, renamed to test_codingConventions, speed up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

disp('Checking all <a href="https://gitlab.lrz.de/cps/cora/-/wikis/home/Coding-conventions">CORA coding conventions</a>..')
assert(aux_checkDocstrings(), 'Some docstrings or file formattings are not following the CORA coding conventions.');

% test completed
res = true;

end


% Auxiliary functions -----------------------------------------------------

function res = aux_checkDocstrings()
% check docstrings and file formatting 
% (everything that requires reading the file content)

coraroot = CORAROOT; % read once for speed
files = [
    findfiles([coraroot filesep 'app']);
    findfiles([coraroot filesep 'contDynamics']);
    findfiles([coraroot filesep 'contSet']);
    findfiles([coraroot filesep 'converter']);
    findfiles([coraroot filesep 'discrDynamics']);
    findfiles([coraroot filesep 'examples']);
    findfiles([coraroot filesep 'global']);
    findfiles([coraroot filesep 'hybridDynamics']);
    findfiles([coraroot filesep 'matrixSet']);
    findfiles([coraroot filesep 'models']);
    findfiles([coraroot filesep 'nn']);
    findfiles([coraroot filesep 'specification']);
    findfiles([coraroot filesep 'unitTests']);
];

% exclude:
% - app/auxiliary files
% files = excludefiles(files, ['app' filesep 'auxiliary']);
% - thirdparty files
files = excludefiles(files, ['global' filesep 'thirdparty']);
% - models: generated/converted/auxiliary functions (see also resetCORA)
files = excludefiles(files, ['models' filesep 'auxiliary']);
files = excludefiles(files, ['models' filesep 'CoraConverted']);
files = excludefiles(files, ['models' filesep 'powerSystemsConverted']);
files = excludefiles(files, ['models' filesep 'spaceExConverted']);
% - repeatability package: main file (no author block)
files = excludefiles(files, ['unitTests' filesep 'ci' filesep 'repeatability-template' filesep 'code'],'main.m');

% iterate through all files
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

        % should have at least one example call
        cnt = 0;
        while lcnt <= length(lines) && startsWith(lines{lcnt},'%') && ~isEmptyCommentLine(lines{lcnt})
            if ~isempty(regexp(lines{lcnt}, sprintf('[\\.= ]%s(\\(.*|)$',filename), 'once','emptymatch'))
                cnt = cnt + 1;
            end
            lcnt = lcnt + 1;
        end
        if cnt == 0
            issues{end+1} = 'At least one example syntax call should be stated.';
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
            issues{end+1} = ['END OF ' 'CODE is written inconsistently.'];
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
        if startsWith(file.folder, [coraroot filesep 'unitTests' filesep])
            if contains(filetext, ['close' 'all'])
                issues{end+1} = "Avoid 'close all' in unit tests.";
            end
        end

        % avoid more than 2 empty lines
        cnt = 0; MAX_EMPTY = 2;
        for lcnt = 1:length(lines)
            line = lines{lcnt};
            if all(line == ' ' | line == '%') % compute directly for speed
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

        % find undocummented code
        cnt = 0; lcnt = 1; MAX_LINES_WO_COMMENT = 25; 
        while lcnt < numel(lines)
            line = lines{lcnt};
            if ~contains(line,'%')
                cnt = cnt + 1;
            else
                if cnt > MAX_LINES_WO_COMMENT
                    issues{end+1} = sprintf('Undocumented code (e.g., <a href="matlab:opentoline(''%s'',%i)">lines %i-%i</a>).',file_path,lcnt-cnt,lcnt-cnt,lcnt-1);
                    break
                end

                % long undocumented code is ok for model files
                if contains(line,'% model') % <- marker for models
                    lcnt = lcnt + 1;
                    % skip to next comment
                    while ~contains(lines{lcnt},'%')
                        lcnt = lcnt + 1;
                    end
                end

                % reset counter
                cnt = 0;
            end

            % increase counter
            lcnt = lcnt + 1;
        end

        % check see also for function_ files

        if endsWith(file.name, '_.m') && contains(file_path,'contSet')
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

        % check error(...) -> CORAerror(...)
        probCall = 'error(';
        allowedCalls = {'CORAerror(','.error(','_error(','yalmiperror('};
        if ~ismember(filename,{'CORAerror','test_codingConventions'}) && ...
                ~aux_checkFunctionCall(filetext,probCall,allowedCalls)
            issues{end+1} = "Please replace error(...) calls with CORAerror(id, ...)";
        end

        % check CORAerror(...) -> throw(CORAerror(...))
        probCall = 'CORAerror(';
        allowedCalls = {'throw(CORAerror(','throwAsCaller(CORAerror('};
        if ~ismember(filename,{'CORAerror','test_codingConventions'}) && ...
                ~aux_checkFunctionCall(filetext,probCall,allowedCalls)
            issues{end+1} = "Please replace CORAerror(...) calls with throw(CORAerror(id, ...))";
        end

        % check warning(...) -> CORAwarning(...)
        probCall = 'warning(';
        allowedCalls = {'CORAwarning(','warning(''on''','warning(''off''','warning()','warning(w)','warning(warOrig)'};
        if ~ismember(filename,{'CORAwarning','test_codingConventions'}) && ...
                ~aux_checkFunctionCall(filetext,probCall,allowedCalls)
            issues{end+1} = "Please replace warning(...) calls with CORAwarning(id, ...)";
        end

        % check linprog(...) -> CORAlinprog(...)
        probCall = 'linprog(';
        allowedCalls = {'CORAlinprog(','intlinprog('};
        if ~ismember(filename,{'CORAlinprog','test_codingConventions'}) && ...
            ~aux_checkFunctionCall(filetext,probCall,allowedCalls)
            issues{end+1} = "Please replace linprog(...) calls with CORAlinprog(problem)";
        end

        % check quadprog(...) -> CORAquadprog(...)
        probCall = 'quadprog(';
        allowedCalls = {'CORAquadprog('};
        if ~ismember(filename,{'CORAquadprog','test_codingConventions'}) && ...
            ~aux_checkFunctionCall(filetext,probCall,allowedCalls)
            issues{end+1} = "Please replace quadprog(...) calls with CORAquadprog(problem)";
        end
    
        % spelling
        if contains(lower(filetext),['un' 'kown'])
            issues{end+1} = "Misspelled 'unknown'.";
        end

        % check evParams -> options.nn
        if ~strcmp(filename,'test_codingConventions') && contains(filetext, 'evParams')
            issues{end+1} = 'With appropriate changes, please replace evParams with options.nn.';
        end

        % test deprecated error messages
        if contains(filetext,['CORA:' 'testFailed'])
            issues{end+1} = "Please replace 'testFailed' CORA errors with assert/assertLoop.";
        end
        if contains(filetext,['CORA:' 'tooManyInputArgs']) || contains(filetext,['CORA:' 'notEnoughInputArgs'])
            issues{end+1} = "Please replace number of input argument checks with narginchk.";
        end

        if contains(filetext,['>>>>' '>>>']) || contains(filetext,['<<<<' '<<<' ' HEAD:'])
            issues{end+1}  = "Unresolved merge conflict in " + file.name;
        end

        % check private function
        if contains(file_path, [filesep 'private' filesep]) && ~startsWith(filename,'priv_')
            issues{end+1} = "Private functions should have the prefix 'priv_'";
        end

        % everything tested -----------------------------------------------

        everythingTested = true;

    end

    % check if everything was tested
    if ~everythingTested
        issues{end+1} = '... and potentially some others.';
    end

    % display information
    all_good = isempty(issues);
    if ~all_good
        fprintf('- <a href="matlab:open %s">%s</a>\n', file_path,file_path)
        disp("  Issues:")
        fprintf("  - %s\n", string(issues'))
        fprintf('  (Maybe try to <a href="matlab:fix_docstring(''%s'')">fix the docstring automatically</a>?)\n', file_path)
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
    fprintf('Found %i file(s) violating the <a href="https://gitlab.lrz.de/cps/cora/-/wikis/home/Coding-conventions">CORA coding conventions</a>.\n', issue_count);
end

end

function res = aux_checkFunctionCall(filetext, probCall,allowedCalls)
% checks if all problematic calls of a function are replaced by the
% allowed calls

res = true;

% find indices with potentially problematic calls
idx = strfind(filetext,probCall);

% check each occurence
for pos = idx
    % check if matches with allowed calls
    resvec = false(1,numel(allowedCalls));
    for i=1:numel(allowedCalls)
        allowedCall = allowedCalls{i};

        % find start position of problematic call in allowed call
        pos_i = strfind(allowedCall,probCall);

        % extract string including enough text before and after to check if
        % it is a allowed call
        potProbString = filetext((pos-pos_i):(pos-pos_i+numel(allowedCall)+1));

        % check if it is a allowed call
        resvec(i) = contains(potProbString,allowedCall);
        if any(resvec)
            break;
        end
    end

    % if none of the allowed calls matched, we found a counterexample
    res = any(resvec);
    if ~res
        break
    end
end

end

% ------------------------------ END OF CODE ------------------------------
