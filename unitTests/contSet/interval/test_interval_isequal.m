function res = test_interval_isequal
% test_interval_isequal - unit test function of isequal
%
% Syntax:  
%    res = test_interval_isequal
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
upper = [4;   2;  6;  3;  8];
Int1 = interval(lower, upper);
upper = [4;   2;  6;  2;  8];
Int2 = interval(lower, upper);
Int3 = Int1;

% check if all points are in interval
res = ~isequal(Int1,Int2) && isequal(Int1,Int3);

if res
    disp('test_isequal successful');
else
    disp('test_isequal failed');
end

%------------- END OF CODE --------------