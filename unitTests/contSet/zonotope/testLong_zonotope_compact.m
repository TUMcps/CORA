function res = testLong_zonotope_compact
% testLong_zonotope_compact - unit test function of compact
%    this encompasses checking the function nonzeroFilter
%
% Syntax:
%    res = testLong_zonotope_compact
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

% Authors:       Matthias Althoff, Mark Wetzlinger
% Written:       26-July-2016
% Last update:   09-August-2020 (MW, enhance randomness)
% Last revision: 29-July-2023 (MW, rename '...compact')

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;
tol = 1e-12;

% number of tests
nrTests = 1000;

% compare randomly generated zonotopes
for test=1:nrTests

    % random dimension, random number of generators
    n = randi(10);
    nrOfGens = randi([2*n,5*n]);
    
    % center, generator matrix without any all-zero generators
    c = zeros(n,1);
    G_nozeros = 1+rand(n,nrOfGens);

    % zonotope without all-zero generators
    Z_nozeros = zonotope(c,G_nozeros);

    % representation without zero-generators
    Z_compact = compact_(Z_nozeros,'zeros',tol);
    
    % since no zero generators, results has to be the same as before
    assert(compareMatrices(Z_compact.c,c) ...
        && compareMatrices(Z_compact.G,G_nozeros));
    
    % append zero generators
    G_withzeros = [G_nozeros, zeros(n,randi([5,25],1,1))];
    % shuffle matrix
    G_withzeros = G_withzeros(:,randperm(size(G_withzeros,2)));
    Z_withzeros = zonotope(c,G_withzeros);
    
    % representation without zero-generators
    Z_compact = compact_(Z_withzeros,'zeros',tol);
    
    % result has to be the same as original zonotope
    assert(isequal(Z_compact,Z_nozeros,tol));
    assert(compareMatrices(Z_compact.c,c) ...
        && compareMatrices(Z_compact.G,G_nozeros));

end

% ------------------------------ END OF CODE ------------------------------
