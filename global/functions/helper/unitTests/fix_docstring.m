function fix_docstring(file_path)
% fix_docstring - this function aims to fix the docstring of a given file
%
% Syntax:
%    res = fix_docstring(file_path)
%
% Inputs:
%    file_path - char
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_docstring

% Authors:       Tobias Ladner
% Written:       17-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

fprintf("Trying to fix file: %s\n", file_path)
disp("Please rerun test_docstring afterwards.")

% read content of file
filetext = fileread(file_path);

% split lines
lines = splitlines(filetext);
changed = false;

for h = 1 % hack to use break within code checks

    % initial checks --------------------------------------------------

    if isempty(lines)
        disp("Empty file. Please provide basic docstring structure.")
        break
    end

    if ~startsWith(lines{1}, 'function') && ~startsWith(filetext, 'classdef')
        disp("Given file is a script file.")
        break
    end

    % init line counter
    lcnt = 2;

    % unable to fix description automatically

    lcnt = lcnt + 2;

    % test syntax block in docstring ----------------------------------

    % iterate to syntax block
    while lcnt <= length(lines) && ~startsWith(lines{lcnt}, '% Syntax')
        lcnt = lcnt + 1;
    end

    if lcnt > length(lines)
        disp("Syntax block is missing. Please provide basic docstring structure.")
        break
    end

    % update syntax line
    lines{lcnt} = '% Syntax:';
    changed = true;

    if ~isEmptyCommentLine(lines{lcnt-1})
        % add empty line prior to syntax block
        lines = [lines(1:lcnt-1); {'%'}; lines(lcnt:end)];
        lcnt = lcnt + 1;
    end

    % unable to add exemplary call

    % test inputs block in docstring ----------------------------------

    % iterate to inputs block
    while lcnt <= length(lines) && ~startsWith(lines{lcnt}, '% Input')
        lcnt = lcnt + 1;
    end

    if lcnt > length(lines)
        disp("Inputs block is missing. Please provide basic docstring structure.")
        break
    end

    % update inputs line
    lines{lcnt} = '% Inputs:';

    if ~isEmptyCommentLine(lines{lcnt-1})
        % add empty line prior to syntax block
        lines = [lines(1:lcnt-1); {'%'}; lines(lcnt:end)];
        lcnt = lcnt + 1;
    end

    % unable to add inputs description

    % test outputs block in docstring ---------------------------------

    % iterate to output block
    while lcnt <= length(lines) && ~startsWith(lines{lcnt}, '% Output')
        lcnt = lcnt + 1;
    end

    if lcnt > length(lines)
        disp("Outputs block is missing. Please provide basic docstring structure.")
        break
    end

    % update inputs line
    lines{lcnt} = '% Outputs:';

    if ~isEmptyCommentLine(lines{lcnt-1})
        % add empty line prior to syntax block
        lines = [lines(1:lcnt-1); {'%'}; lines(lcnt:end)];
        lcnt = lcnt + 1;
    end

    % unable to add inputs description

    % test authors block in docstring ---------------------------------

    % iterate to author block
    while lcnt <= length(lines) && ~startsWith(lines{lcnt}, '% Author')
        lcnt = lcnt + 1;
    end

    if lcnt > length(lines)
        disp('Author block is missing.');
        break
    end

    % check if previous line is empty
    if ~isempty(lines{lcnt-1})
        lines = [lines(1:lcnt-1); {''}; lines(lcnt:end)];
        lcnt = lcnt + 1;
    end

    % fix authors text
    lines{lcnt} = regexprep(lines{lcnt},'% Author(s|): *','% Authors:       ');
    lcnt = lcnt + 1;

    % fix written text
    lines{lcnt} = regexprep(lines{lcnt},'% Written: *','% Written:       ');
    lines{lcnt} = aux_fixAuthorDateLine(lines{lcnt});
    lcnt = lcnt + 1;

    % fix last update text
    lines{lcnt} = regexprep(lines{lcnt},'% Last update: *','% Last update:   ');

    while lcnt <= length(lines) && ~startsWith(lines{lcnt},'% Last revision:')
        lines{lcnt} = aux_fixAuthorDateLine(lines{lcnt});
        lcnt = lcnt + 1;
    end

    % fix last revision text
    if lcnt > length(lines)
        disp('Last revision block is missing.')
        break
    end


    lines{lcnt} = regexprep(lines{lcnt},'% Last revision: *','% Last revision: ');

    while lcnt <= length(lines) && ~isEmptyLine(lines{lcnt})
        lines{lcnt} = aux_fixAuthorDateLine(lines{lcnt});
        lcnt = lcnt + 1;
    end

    % test BEGIN CODE -------------------------------------------------

    % iterate to BEGIN CODE
    while lcnt <= length(lines) && ~contains(lines{lcnt}, 'BEGIN CODE')
        lcnt = lcnt + 1;
    end

    if lcnt > length(lines)
        break
    end

    lines{lcnt} = '% ------------------------------ BEGIN CODE -------------------------------';

    % check if previous line is empty
    if ~isEmptyLine(lines{lcnt-1})
        % add empty line
        lines = [lines(1:lcnt-1); {''}; lines(lcnt:end)];
        lcnt = lcnt + 1;
    end

    % only single empty line
    while isEmptyLine(lines{lcnt-2})
        % remove empty line
        lines = [lines(1:lcnt-3); lines(lcnt-1:end)];
        lcnt = lcnt + 1;
    end

    % check if subsequent line is empty
    if ~isEmptyLine(lines{lcnt+1})
        % add empty line
        lines = [lines(1:lcnt); {''}; lines(lcnt+1:end)];
        lcnt = lcnt + 1;
    end

    % test END CODE -------------------------------------------------

    % iterate to END CODE
    while lcnt <= length(lines) && ~contains(lines{lcnt}, ['END OF ', 'CODE'])
        lcnt = lcnt + 1;
    end

    if lcnt > length(lines)
        break
    end

    lines{lcnt} = ['% ------------------------------ END OF ', 'CODE ------------------------------'];

    % check if previous line is empty
    if ~isEmptyLine(lines{lcnt-1})
        % add empty line
        lines = [lines(1:lcnt-1); {''}; lines(lcnt:end)];
        lcnt = lcnt + 1;
    end

    % check if subsequent line is empty
    if lcnt + 1 > length(lines) || isEmptyLine(lines{lcnt+1})
        % add empty line
        lines = [lines(1:lcnt); {''}; lines(lcnt+1:end)];
        lcnt = lcnt + 1;
    end

    % additional checks -----------------------------------------------

    % auxiliary functions
    lcnt = 2;
    while lcnt <= length(lines) && ~startsWith(lines{lcnt}, 'function')
        lcnt = lcnt + 1;
    end

    if lcnt <= length(lines)
        % first auxiliary function found

        % iterate back to find '% Auxiliary functions' block
        lcntp = lcnt;
        while lcntp > 0 && ~startsWith(lower(lines{lcntp}), '% aux')
            lcntp = lcntp - 1;
        end

        if lcntp == 0
            fprintf("Unable to find '% Auxiliary functions' block. Maybe there is a nested function at line %d? \n", lcnt)
        else
            lines{lcntp} = '% Auxiliary functions -----------------------------------------------------';

            % check if previous two lines are empty
            if ~isEmptyLine(lines{lcntp-1})
                % add empty line
                lines = [lines(1:lcntp-1); {''}; lines(lcntp:end)];
                lcnt = lcnt + 1;
                lcntp = lcntp + 1;
            end
            if ~isEmptyLine(lines{lcntp-2})
                % add empty line
                lines = [lines(1:lcntp-2); {''}; lines(lcntp-1:end)];
                lcnt = lcnt + 1;
                lcntp = lcntp + 1;
            end

            % check if subsequent line is empty
            if ~isEmptyLine(lines{lcntp+1})
                % add empty line
                lines = [lines(1:lcntp); {''}; lines(lcntp+1:end)];
                lcnt = lcnt + 1;
            end
        end

        % auxiliary functions should start with 'aux_'
        while lcnt <= length(lines)
            if startsWith(lines{lcnt}, 'function ') && ~contains(lines{lcnt}, 'aux_')
                % replace all functions calls with 'aux_' prefix

                % gather function name
                fname = lines{lcnt};
                fname = fname(10:end); % remove 'function '
                
                % remove paranthesis
                idxP = strfind(fname, '(');
                while isempty(idxP) || lcnt == length(lines)
                    lcnt = lcnt + 1;
                    fname = lines{lcnt};
                    idxP = strfind(fname, '(');
                end
                if ~isempty(idxP)
                    fname = fname(1:idxP-1); % remove after '('
                else
                    fprintf("Unable to find function name: '%s'\n", lines{lcnt});
                    continue
                end

                % remove everything up to equal sign
                idxEq = strfind(fname, '=');
                if ~isempty(idxEq)
                    fname = fname(idxEq+1:end); % remove up to '='
                end
                
                fname = strtrim(fname); % remove spaces
                if contains(fname, ' ')
                    fprintf("Unable to find function name: '%s'\n", lines{lcnt});
                end

                if startsWith(fname,'aux_')
                    % all good
                    continue
                end

                aux_fname = ['aux_', fname];

                for lcntf = 1:length(lines)
                    if contains(lines{lcntf}, [fname '(']) || contains(lines{lcntf}, ['@' fname])
                        lines{lcntf} = strrep(lines{lcntf}, fname, aux_fname);
                    end
                end

            end
            lcnt = lcnt + 1;
        end
    end

    % fix more than two empty lines
    lcnt = 1; cnt = 0; MAX_EMPTY = 2;
    while lcnt <= length(lines)
        if isEmptyCommentLine(lines{lcnt}) || isEmptyLine(lines{lcnt})
            cnt = cnt + 1;
        else
            if cnt > MAX_EMPTY
                lines = [lines(1:(lcnt-cnt+MAX_EMPTY-1)); lines(lcnt:end)];
                lcnt = lcnt - cnt + MAX_EMPTY;
            end

            % reset counter
            cnt = 0;
        end

        lcnt = lcnt + 1;
    end

end

% remove empty lines at the end (gets implicitly added in writecell)
lcnt = length(lines);
while isEmptyLine(lines{lcnt})
    lines = lines(1:end-1);
    lcnt = lcnt - 1;
end

% write lines back to file
writecell(lines, file_path, 'FileType', 'text', 'QuoteStrings', 'none');

end


% Auxiliary functions -----------------------------------------------------

function line = aux_fixAuthorDateLine(line)

    if length(line) < 17
        % unable to fix
        return
    end

    % fix empty space
    if line(17) ~= ' '
        line = [line(1:16) ' ' line(17:end)];
    end

    if length(line) < 19
        line = [line(1:17) '---'];
        return
    end

    % date starts with single digit
    if ~isempty(regexp(line(18:19),'[0-9]-', 'once'))
        line = [line(1:17) '0' line(18:end)];
    end

    % fix short month
    line = strrep(line,'-Jan-','-January-');
    line = strrep(line,'-Feb-','-February-');
    line = strrep(line,'-Mar-','-March-');
    line = strrep(line,'-Apr-','-April-');
    % line = strrep(line,'-May-','-May-');
    line = strrep(line,'-Jun-','-June-');
    line = strrep(line,'-Jul-','-July-');
    line = strrep(line,'-Aug-','-August-');
    line = strrep(line,'-Sep-','-September-');
    line = strrep(line,'-Oct-','-October-');
    line = strrep(line,'-Nov-','-November-');
    line = strrep(line,'-Dec-','-December-');

    % fix comment
    idx = regexp(line,'([A-Z][A-Z]:','once');
    if ~isempty(idx)
        line(idx+3) = ',';
    end

end

% ------------------------------ END OF CODE ------------------------------
