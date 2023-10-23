function res = testLong_zonotope_convHull
% testLong_zonotope_convHull - unit test function of convHull
%
% Syntax:
%    res = testLong_zonotope_convHull
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

% assume true
res = true;

% number of tests
nrTests = 25;
% number of sampled points
nrRandPoints = 100;

% box has to be the same as conversion to interval
for i=1:nrTests

    % random dimension
    n = randi(8);
    % number of generators (keep low to speed up containment check)
    nrGens = randi([n,2*n]);

    % init two random zonotopes
    Z1 = zonotope.generateRandom('Dimension',n,'NrGenerators',nrGens);
    Z2 = zonotope.generateRandom('Dimension',n,'NrGenerators',nrGens);

    % compute convex hull
    Z = convHull(Z1,Z2);

    % sample points from individual capsules
    p1 = randPoint(Z1,nrRandPoints);
    p2 = randPoint(Z1,nrRandPoints);

    % compute linear combination of pairs
    lambda = rand(1,nrRandPoints);
    p = lambda .* p1 + (1-lambda) .* p2;

    % check if those points are contained in the convex hull
    if ~all(contains(Z,[p1 p2 p]))
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
