function res = testLong_polytope_generateRandom
% testLong_polytope_generateRandom - random tests for the
%    generation of randomized polytopes
%
% Syntax:
%    res = testLong_polytope_generateRandom()
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       01-December-2022
% Last update:   12-December-2022 (MW, add emptiness check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 100;

for i = 1:nrTests

    % random dimension
    n = randi(6);

    % random number of constraints
    nrCon = randi([ceil(n/2),2*n]);

    % create random polytopes

    % only dimension
    P = polytope.generateRandom('Dimension',n);
    assertLoop(dim(P) == n,i)

    % should not be empty
    assertLoop(~representsa(P, 'emptySet'),i)

    % dimension and number of constraints
    P = polytope.generateRandom('Dimension',n,'NrConstraints',nrCon);
    assertLoop(dim(P) == n,i)
    assertLoop(size(P.A,1) == nrCon,i)
    assertLoop(~representsa(P, 'emptySet'),i)

    % not enough constraints for bounded polytope
    P = polytope.generateRandom('Dimension',n,'NrConstraints',n);
    assertLoop(~isBounded(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    % boundedness
    P = polytope.generateRandom('IsBounded',false);
    assertLoop(~isBounded(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    P = polytope.generateRandom('Dimension',n,'IsBounded',false);
    assertLoop(~isBounded(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    P = polytope.generateRandom('Dimension',n,'NrConstraints',2*n,...
        'IsBounded',false);
    assertLoop(~isBounded(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    % degeneracy
    P = polytope.generateRandom('IsDegenerate',true);
    assertLoop(~isFullDim(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    P = polytope.generateRandom('Dimension',n,'IsDegenerate',true);
    assertLoop(~isFullDim(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    P = polytope.generateRandom('Dimension',n,'NrConstraints',2*n,'IsDegenerate',true);
    assertLoop(~isFullDim(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)

    % polytope with 2 constraints can be unbounded, but degenerate
    P = polytope.generateRandom('Dimension',n,'NrConstraints',2,'IsDegenerate',true);
    assertLoop(~isFullDim(P),i)
    assertLoop(~representsa(P, 'emptySet'),i)
end

% ------------------------------ END OF CODE ------------------------------
