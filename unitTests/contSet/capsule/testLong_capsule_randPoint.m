function res = testLong_capsule_randPoint
% testLong_capsule_randPoint - unit test function of randPoint
%
% Syntax:
%    res = testLong_capsule_randPoint
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: @capsule > containsPoint.m
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       17-September-2019
% Last update:   12-March-2021 (MW, more thorough testing)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

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
    res = all(contains(C,p));
    
    if ~res
        break;
    end
end

% ------------------------------ END OF CODE ------------------------------
