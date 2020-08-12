function res = test_capsule_randPoint
% test_capsule_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = test_capsule_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: @capsule > containsPoint.m
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      17-Sep-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% instantiate capsule
C = capsule([1;1],[0;1],0.5);

% generate random points
numPoints = 100;
p = zeros(dim(C),numPoints);
for i=1:numPoints
    p(:,i) = randPoint(C);
end

% check if all random points inside capsule
res = all(containsPoint(C,p));

if res
    disp('test_randPoint successful');
else
    disp('test_randPoint failed');
end

%------------- END OF CODE --------------