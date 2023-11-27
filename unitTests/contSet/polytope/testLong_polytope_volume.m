function res = testLong_polytope_volume
% testLong_polytope_volume - unit test function of volume
%
% Syntax:
%    res = testLong_polytope_volume
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
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% generate zonotopes in low dimensions, compute halfspace representation
% and compare values of the two volume functions: since the conversion from
% zonotopes to polytopes is exact, both volumes should be identical up to
% numerical precision

nrTests = 25;

for i=1:nrTests
    
    % random dimension
    n = randi([2,4]);

    % random number of generators
    nrGens = randi([n,2*n]);

    % generate random zonotope, where the vertex enumeration of the
    % converted polytope is successful
    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',nrGens);
    % scale to unit box to avoid numerical imprecisions
    Z = enlarge(Z,1./rad(interval(Z)));
    P = polytope(Z);
    try
        vertices(P);
    catch ME
        continue
    end

    % compute volumes
    val_Z = volume(Z);
    val_P = volume(P);

    % compare values
    if ~withinTol(val_Z,val_P,1e-5)
        throw(CORAerror('CORA:testFailed'));
    end

end


% ------------------------------ END OF CODE ------------------------------
