function res = testLong_updateCORApath
% testLong_updateCORApath - tests updating the CORA path
%
% Syntax:
%    res = testLong_updateCORApath()
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
% See also: none

% Authors:       Tobias Ladner
% Written:       23-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% turn off warning while re-setting path
w = warning;
warning off;
    
% test run
updateCORApath();

% remove directory from path and check if is on path
paths = {
    [CORAROOT filesep 'contSet'];
    [CORAROOT filesep 'contDynamics'];
    [CORAROOT filesep 'global' filesep 'classes'];
};

for i=1:size(paths,1)
    % remove directory from matlab path
    path = paths{i};
    rmpath(path);

    % restore
    updateCORApath();

    % check if on matlab path
    pathCell = regexp(path, pathsep, 'split');
    assert(any(strcmp(path, pathCell)));
end


% restore warning
warning(w)

% gather results
res = true;

% ------------------------------ END OF CODE ------------------------------
