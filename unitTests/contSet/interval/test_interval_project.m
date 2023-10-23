function res = test_interval_project
% test_interval_project - unit test function of project
%
% Syntax:
%    res = test_interval_project
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

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% create interval
lower = [-3; -9; -4; -7; -1];
upper = [4; 2; 6; 3; 8];
Int = interval(lower, upper);

% project interval
dimensions = [1;3;4];
I_proj1 = project(Int,dimensions);

% logical indexing
dimensions = [true false true true false];
I_proj2 = project(Int,dimensions);

% true solution
lower = [-3; -4; -7];
upper = [ 4;  6;  3];
I_true = interval(lower, upper);

% check if all points are in interval
res = isequal(I_proj1,I_true) && isequal(I_proj2,I_true);

% ------------------------------ END OF CODE ------------------------------
