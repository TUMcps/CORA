function res = testLong_capsule_convHull
% testLong_capsule_convHull - unit test function of convHull
%
% Syntax:
%    res = testLong_capsule_convHull
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
% Written:       23-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nrTests = 1000;
nrRandPoints = 100;

for i=1:nrTests
    % random dimension
    n = randi([2,8]);

    % init two random capsules
    C1 = capsule.generateRandom('Dimension',n);
    C2 = capsule.generateRandom('Dimension',n);

    % compute convex hull
    C = convHull(C1,C2);

    % sample points from individual capsules
    p1 = randPoint(C1,nrRandPoints);
    p2 = randPoint(C1,nrRandPoints);

    % compute linear combination of pairs
    lambda = rand(1,nrRandPoints);
    p = lambda .* p1 + (1-lambda) .* p2;

    % check if those points are contained in the convex hull
    if ~all(contains(C,[p1 p2 p]))
        res = false; break
    end

end

% ------------------------------ END OF CODE ------------------------------
