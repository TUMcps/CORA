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
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% turn off warning while re-setting path
w = warning;
warning off;

% add all files from corapath to path
corapath = CORAROOT;
rmpath(genpath(corapath))
addpath(genpath(corapath))
rehash path
% savepath

% restore warning
warning(w)

% restore class definitions (e.g. after switching git branch)
clear classes

end

% ------------------------------ END OF CODE ------------------------------
