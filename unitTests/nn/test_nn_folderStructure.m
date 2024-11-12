function res = test_nn_folderStructure()
% test_nn_folderStructure - tests the folder structure of all nn-related files
%    as nn stuff was moved from <CORAROOT>/global/classes/nn to <CORAROOT>/nn
%
% Syntax:
%    res = test_nn_folderStructure()
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
% See also: -

% Authors:       Tobias Ladner
% Written:       22-July-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% find all directories
files = dir([CORAROOT filesep '**' filesep]);
dirs = files([files.isdir]);

% specify allowed nn paths
allowedNNPaths = {
    [CORAROOT filesep 'nn'];
    [CORAROOT filesep 'examples' filesep 'nn'];
    [CORAROOT filesep 'examples' filesep 'website' filesep 'nn'];
    [CORAROOT filesep 'global' filesep 'functions' filesep 'helper' filesep 'nn'];
    [CORAROOT filesep 'models' filesep 'Cora' filesep 'nn'];
    [CORAROOT filesep 'unitTests' filesep 'nn'];
    [CORAROOT filesep 'manual']; % don't care about manual folder structure
};

% find faulty paths
idxFaulty = arrayfun(@(path) ...
    ( ...
        contains([path.folder filesep  path.name], [filesep 'nn' filesep]) || ...
        contains([path.folder filesep  path.name], [filesep 'neuralNetworks' filesep]) ...
    ) && ...
    ~contains([path.folder filesep  path.name], allowedNNPaths), ...
    dirs ...
);
dirsFaulty = dirs(idxFaulty);

% check if no faulty directories exist
if ~isempty(dirsFaulty)
    % list faulty dirs
    disp('Faulty directories:')
    for i=1:numel(dirsFaulty)
        fprintf('- %s\n', [dirsFaulty(i).folder filesep dirsFaulty(i).name])
    end
end
assert(isempty(dirsFaulty));

end

% ------------------------------ END OF CODE ------------------------------
