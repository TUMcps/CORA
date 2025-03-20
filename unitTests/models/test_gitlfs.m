function res = test_gitlfs()
% test_gitlfs - tests if all files were properly downloaded using git lfs
%
% Syntax:
%    res = test_gitlfs()
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
% Written:       14-February-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% load .mat file (synched via git lfs)
model = "gitlfs.mat";
load(model,'gitlfs')

res = gitlfs;

% ------------------------------ END OF CODE ------------------------------
