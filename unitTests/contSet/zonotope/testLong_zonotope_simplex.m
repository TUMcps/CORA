function res = testLong_zonotope_simplex
% testLong_zonotope_simplex - unit test function of simplex
%
% Syntax:
%    res = testLong_zonotope_simplex
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
nrTests = 100;

% box has to be the same as conversion to interval
for i=1:nrTests

    % random dimension
    n = randi(10);
    % number of generators
    nrGens = randi(5*n);

    % init two random zonotopes
    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',nrGens);

    % convert to simplex
    P = simplex(Z);

    % check if zonotope is contained in simplex
    if ~contains(P,Z,'exact',1e-12)
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
