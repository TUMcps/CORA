function resetCORA()
% resetCORA - resets CORA by removing all converted and auxiliary files
%
% Syntax:
%    resetCORA()
%
% Inputs:
%    -
%
% Outputs:
%    -
%

% Authors:       Tobias Ladner
% Written:       12-September-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% delete converted models -------------------------------------------------

aux_deleteFolder([CORAROOT filesep 'models' filesep 'auxiliary']);
aux_deleteFolder([CORAROOT filesep 'models' filesep 'CoraConverted']);
aux_deleteFolder([CORAROOT filesep 'models' filesep 'powerSystemsConverted']);
aux_deleteFolder([CORAROOT filesep 'models' filesep 'SpaceExConverted']);

% reset macros ------------------------------------------------------------

aux_resetMacro('CHECKS_ENABLED')
aux_resetMacro('CORA_WARNINGS_ENABLED')

% CORA reset --------------------------------------------------------------

end


% Auxiliary functions -----------------------------------------------------

function aux_deleteFolder(path)

    if exist(path,'dir')
        % temporarily disable warnings
        % rmdir warns user if a path was on the MATLAB path
        w = warning();
        warning off;
        % remove folders and subfolders
        res = rmdir(path, 's');
        warning(w); % restore warning
    
        if ~res
            CORAwarning('CORA:global','Unable to remove path: %s', path);
        end
    end

end

function aux_resetMacro(MACRONAME)

% read file
file_path = [CORAROOT filesep 'global' filesep 'macros' filesep MACRONAME '.m'];
filetext = fileread(file_path);
lines = splitlines(filetext);

% find line
for lcnt = 1:length(lines)
    if contains(lines{lcnt},'res = ')
        lines{lcnt} = strrep(lines{lcnt}, 'false', 'true');
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

% update path
updateCORApath();

% check if successful
if ~eval(MACRONAME)
    CORAwarning('CORA:global','Unable to reset CHECKS_ENABLED.')
end

end

% ------------------------------ END OF CODE ------------------------------
