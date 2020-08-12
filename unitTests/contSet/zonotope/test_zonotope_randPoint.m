function res = test_zonotope_randPoint
% test_zonotope_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = test_zonotope_randPoint
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

% create zonotope
Z = zonotope([0, 2,-1;
              1, 1, 1]);

% obtain zonotope without zeros
numPoints = 100;
p = zeros(dim(Z),numPoints);
for i=1:numPoints
    p(:,i) = randPoint(Z);
end

% check result
res = all(containsPoint(Z,p));


if res
    disp('test_randPoint successful');
else
    disp('test_randPoint failed');
end

%------------- END OF CODE --------------
