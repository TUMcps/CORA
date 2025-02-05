function res = updateUnitTests()
% updateUnitTests - this function updates the unit tests from 
%     'res(end+1) = <cond>' to 'assert(<cond>)'
%
% Syntax:
%    res = updateUnitTests()
%
% Inputs:
%    -
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: test_codingConventions

% Authors:       Tobias Ladner
% Written:       28-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

files = findfiles([CORAROOT filesep 'unitTests'],true,'test');

fprintf('Processing %i files:\n', numel(files))

for f = 1:numel(files)

    % get file
    file = files(f);
    file_path = [file.folder filesep file.name];
    fprintf('- %i/%i: <a href="matlab:open %s">%s</a>\n',f,numel(files),file_path,file.name(1:end-2))

    % skip files in ./unitTests
    if contains(file_path,[CORAROOT filesep 'unitTests'])
        continue
    end

    % read content of file
    filetext = fileread(file_path);
    
    % split lines
    lines = splitlines(filetext);

    % iterate through lines and replace as needed
    for lcnt = 1:numel(lines)

        % get line
        line = lines{lcnt};

        % hack to use break within
        for h = 1

            % ignore comments
            if startsWith(strtrim(line),'%')
                break
            end
    
            % init ---
            if contains(line,'resvec = []') || contains(line,'resvec = true(0)')
                line = 'resvec = true;';
                break
            end
            if contains(line,'res = []') || contains(line,'res = true(0)')
                line = 'res = true;';
                break
            end
    
            % assertion ---

            assignments = {'resvec(','res(end+1','res = false','resPartial(','resVec('};
            if any(cellfun(@(assign) contains(line,assign),assignments))
                line = aux_replaceWithAssert(line,lcnt);
                break
            end
            % exclude first line and 'res = true' as it appears at the end
            if lcnt > 1 && ~contains(line,'res = true') && ~contains(line,'res = all') && contains(line,'res = ')
                line = aux_replaceWithAssert(line,lcnt);
                break
            end
            % fix common 

            % manual adaptations ---

            weirdAssignments = {
                'res1','res(1)','resvec(1)','resPartial(1)', ...
                'resi','res(i)','resvec(i)','resPartial(i)'
            };
            if any(cellfun(@(assign) contains(line,assign),weirdAssignments))
                aux_throwManualAdaptationMessage(line,lcnt);
            end
    
    
            % output ---
            if contains(line,'res = all(res)')
                line = 'res = true;';
                break
            end
            if contains(line,'res = all(resvec)')
                line = 'res = true;';
                break
            end

        end

        % fix common miss-corrections
        commonMisCorrection = {'false; break','false;break','false; return','false;return'};
        for i = 1:numel(commonMisCorrection)
            if contains(line,['assert(' commonMisCorrection{i}])
                line = strrep(line,commonMisCorrection{i},'false');
                line = strrep(line,';)',')');
            end
        end

        % save line
        lines{lcnt} = line;

    end

    % check assignment of 'res' after begin of code
    idxBeginCode = find(cellfun(@(line) contains(line,'BEGIN CODE'), lines),1);
    if ~any(cellfun(@(line) contains(line,'res = true'), lines(idxBeginCode:end)))
        % write after begin of code
        lines = [
            lines(1:idxBeginCode);
            {' '};
            {'% assume true'};
            {'res = true;'};
            lines(idxBeginCode+1:end);
        ];
    end
    
    
    % remove empty lines at the end (gets implicitly added in writecell)
    lcnt = numel(lines);
    while isEmptyLine(lines{lcnt})
        lines = lines(1:end-1);
        lcnt = lcnt - 1;
    end
    
    % write lines back to file
    writecell(lines, file_path, 'FileType', 'text', 'QuoteStrings', 'none');

end

disp('Updated all unit tests.')
res = true;

end


% Auxiliary functions -----------------------------------------------------

function line = aux_replaceWithAssert(line,lcnt)
    % assuming line = '    res<variable> = <condition>;'

    % gather relevant indices
    idxEqual = strfind(line,'=');
    idxSemicolon = strfind(line,';');
    idxFirstLetter = strfind(line,'res');
    
    % check if its really 'res =' and not 'S_res', [res,S] or similar
    if (idxFirstLetter(1) > 1 && line(idxFirstLetter(1)-1) == '[') ...
        || contains(line,'...')
        aux_throwManualAdaptationMessage(line,lcnt);
    elseif (idxFirstLetter(1) > 1 && line(idxFirstLetter(1)-1) ~= ' ') ...
        || contains(line,'function res')
        % keep line as is
        % line = line;
    else
        % build new line
        try
            line = sprintf( ...
                '%sassert(%s);%s', ...
                line(1:idxFirstLetter(1)-1), ... % indent
                strtrim(line(idxEqual(1)+1:idxSemicolon(end)-1)), ... % <condition>
                strtrim(line(idxSemicolon(end)+1:end)) ... % comments etc.
            );
        catch ME
            % check manually
            aux_throwManualAdaptationMessage(line,lcnt);
        end
    end
end

function aux_throwManualAdaptationMessage(line,lcnt)
    throw(CORAerror('CORA:specialError',sprintf('Please check line %i manually:\n''%s''',lcnt,line)))
end

% ------------------------------ END OF CODE ------------------------------
