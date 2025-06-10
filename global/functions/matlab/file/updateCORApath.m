function updateCORApath()
% updateCORApath - restores all folders of CORA on the matlab path
%
% Syntax:
%    updateCORApath()
%
% Inputs:
%    -
%
% Outputs:
%    -
%

% Authors:       Tobias Ladner
% Written:       25-November-2022
% Last update:   22-May-2023
%                24-March-2025 (TL, no longer include repeatability template)
%                14-April-2025 (TL, speed up)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% turn off warning while re-setting path
w = warning;
warning off;

% add all files from corapath to path
corapath = CORAROOT;
genedpath = genpath(corapath);
rmpath(genedpath)
addpath(genedpath)

% remove repeatability template to avoid conflicts with other main.m files
rmpath(genpath([corapath filesep 'unitTests' filesep 'ci' filesep 'repeatability-template']))

% rehash path
rehash path

% restore warning
warning(w)

% restore class definitions (e.g., after switching git branch)
clear classes

end

% ------------------------------ END OF CODE ------------------------------
