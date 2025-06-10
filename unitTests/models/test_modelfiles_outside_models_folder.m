function res = test_modelfiles_outside_models_folder()
% test_modelfiles_outside_models_folder - tests if there are any model
%    files stored outside ./models (bad practice as difficult to remove
%    unnecessary files when uploading a repeatability package with size
%    limit; usually, most model files can be deleted, so just removing 
%    ./models is much easier)
%
% Syntax:
%    res = test_modelfiles_outside_models_folder()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: https://gitlab.lrz.de/cps/cora/-/issues/56

% Authors:       Tobias Ladner
% Written:       16-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% settings
modelextensions = {'mat'};
coraroot = CORAROOT; % for speed
modelfolder = [CORAROOT filesep 'models'];

disp('Checking for model files outside ./model folder..');
foundfiles = 0;
for i=1:numel(modelextensions)
    % search files
    files = findfiles(coraroot,true,'',modelextensions{i});

    % filter files ---
    % ... already in ./models
    files = files(arrayfun(@(file) ~contains(file.folder,modelfolder), files));
    % ./unitTests/unitTestsStatus.mat
    files = files(arrayfun(@(file) ~contains(file.name,'unitTestsStatus.mat'), files));

    % any files found?
    foundfiles_i = numel(files);
    foundfiles = foundfiles + foundfiles_i;

    % display found files
    for j=1:foundfiles_i
        file = files(j);
        fprintf('- %s/%s\n',file.folder,file.name)
    end
end
fprintf('Done. Found %i files.\n', foundfiles);

res = foundfiles == 0;

end

% ------------------------------ END OF CODE ------------------------------
