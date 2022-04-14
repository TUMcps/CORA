function res = testLongDuration_capsule_randPoint
% testLongDuration_capsule_randPoint - unit test function of randPoint
%
% Syntax:  
%    res = testLongDuration_capsule_randPoint
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
% Last update:  12-March-2021 (MW, more thorough testing)
% Last revision:---

%------------- BEGIN CODE --------------

numPoints = 100;

% random tests
for n=1:2:20
    
    % instantiate capsule
    C = capsule(randn(n,1),randn(n,1),rand(1));

    % generate random points
    p = zeros(dim(C),numPoints);
    for i=1:numPoints
        p(:,i) = randPoint(C);
    end

    % check if all random points inside capsule
    res = in(C,p);
    
    if ~res
        break;
    end
end


if res
    disp('test_randPoint successful');
else
    disp('test_randPoint failed');
end

%------------- END OF CODE --------------