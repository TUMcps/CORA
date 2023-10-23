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
    if dim(P) ~= n
        res = false; return
    end
    % should not be empty
    if representsa(P, 'emptySet')
        res = false; return
    end

    % dimension and number of constraints
    P = polytope.generateRandom('Dimension',n,'NrConstraints',nrCon);
    if dim(P) ~= n || size(P.A,1) ~= nrCon || representsa(P, 'emptySet')
        res = false; return
    end

    % not enough constraints for bounded polytope
    P = polytope.generateRandom('Dimension',n,'NrConstraints',n);
    if isBounded(P) || representsa(P, 'emptySet')
        res = false; return
    end

    % boundedness
    P = polytope.generateRandom('IsBounded',false);
    if isBounded(P) || representsa(P, 'emptySet')
        res = false; return
    end
    P = polytope.generateRandom('Dimension',n,'IsBounded',false);
    if isBounded(P) || representsa(P, 'emptySet')
        res = false; return
    end
    P = polytope.generateRandom('Dimension',n,'NrConstraints',2*n,...
        'IsBounded',false);
    if isBounded(P) || representsa(P, 'emptySet')
        res = false; return
    end

    % degeneracy
    P = polytope.generateRandom('IsDegenerate',true);
    if isFullDim(P) || representsa(P, 'emptySet')
        res = false; return
    end
    P = polytope.generateRandom('Dimension',n,'IsDegenerate',true);
    if isFullDim(P) || representsa(P, 'emptySet')
        res = false; return
    end
    P = polytope.generateRandom('Dimension',n,'NrConstraints',2*n,...
        'IsDegenerate',true);
    if isFullDim(P) || representsa(P, 'emptySet')
        res = false; return
    end
    % polytope with 2 constraints can be unbounded, but degenerate
    P = polytope.generateRandom('Dimension',n,'NrConstraints',2,...
        'IsDegenerate',true);
    if isFullDim(P) || representsa(P, 'emptySet')
        res = false; return
    end
end

% ------------------------------ END OF CODE ------------------------------
