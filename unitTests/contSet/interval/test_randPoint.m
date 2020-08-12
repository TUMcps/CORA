function res = test_randPoint
% test_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = test_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: @interval > containsPoint.m
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

% compute random points
p = zeros(5,100);
for i=1:100
    p(:,i) = randPoint(Int);
end

% check if all points are in interval
res = all(containsPoint(Int,p));

if res
    disp('test_randPoint successful');
else
    disp('test_randPoint failed');
end

%------------- END OF CODE --------------