function exportAllScriptsForWebsite()
% exportAllScriptsForWebsite - exports the given script to html
%
% Syntax:
%    exportAllScriptsForWebsite()
%
% Inputs:
%    filepath - char, file path
%
% Outputs:
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       30-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% get all files
files = findfiles([CORAROOT '/examples/website'],true,'website_','mlx');

% go through all files
fprintf('Exporting %i files:\n', numel(files))
maxFileName = max(arrayfun(@(file) length(file.name)-4, files));
for i=1:numel(files)
    try
        file = files(i);
        fprintf(['.. %-' num2str(maxFileName) 's : '],file.name(1:end-4));

        % gather path
        filepath = [file.folder filesep file.name];

        % save variables of this function
        vars = who;

        % run file
        matlab.internal.liveeditor.executeAndSave(filepath);

        % clear variables from live script
        clearvars('-except',vars{:})

        % export to html
        exportScriptToHTML(filepath);

        fprintf('done\n')
        
    catch ME
        fprintf('failed (%s)\n', ME.message)
    end
end

disp('Done! Make sure to double-check all exported html files before uploading.')


end

% ------------------------------ END OF CODE ------------------------------
