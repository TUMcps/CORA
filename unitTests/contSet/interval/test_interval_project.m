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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create interval
lower = [-3; -9; -4; -7; -1];
upper = [4; 2; 6; 3; 8];
Int = interval(lower, upper);

% project interval
dimensions = [1;3;4];
Int_project = project(Int,dimensions);

% true solution
lower = [-3; -4; -7];
upper = [ 4;  6;  3];
Int_true = interval(lower, upper);

% check if all points are in interval
res = isequal(Int_project,Int_true);

if res
    disp('test_project successful');
else
    disp('test_project failed');
end

%------------- END OF CODE --------------