function res = testLong_zonotope_isFullDim
% testLong_zonotope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_zonotope_isFullDim
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

% Authors:       Mark Wetzlinger, Adrian Kulmburg
% Written:       12-March-2021
% Last update:   04-February-2024 (AK, added degenerate cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 1000;

for i=1:nrTests

    % random dimension
    n = randi([2,50]);

    % random center
    c = randn(n,1);

    % random number of generators
    gamma = randi([n,100]);

    % random full-rank generator matrix
    G = randn(n,gamma);
    while rank(G) < n
        G = randn(n,gamma);
    end
    
    % instantiate zonotope
    Z = zonotope(c,G);
    
    % check if full-dimensional
    [res_subspace, subspace] = isFullDim(Z);
    assert(isFullDim(Z));
    assert(res_subspace);
    assert((size(subspace,2) == dim(Z)))
end

for i=1:(nrTests/10)
    % random dimension
    n = randi([2,5]);

    % random center
    c = randn(n,1);

    % random number of generators
    gamma = randi([n,10]);

    % random full-rank generator matrix
    G = randn(n,gamma);
    while rank(G) < n
        G = randn(n,gamma);
    end

    % random number of 'degeneracy'
    degeneracy = randi([1,n]);
    [U, Sigma, V] = svd(G);
    for j=1:degeneracy
        Sigma(n-j+1,n-j+1) = 0;
    end
    Gdeg = U * Sigma * V';
    
    % instantiate zonotope
    Z = zonotope(c,Gdeg);
    
    % check if full-dimensional
    res_single = isFullDim(Z);
    [res_subspace, subspace] = isFullDim(Z);

    % transform into polytope, to check if both methods yield the same
    P = polytope(Z);
    [~,P_subspace] = isFullDim(P);

    same_subspace = rank([subspace P_subspace], 1e-6) == rank(subspace);
    assert(~res_single);
    assert(~res_subspace);
    assert(same_subspace);
end

% ------------------------------ END OF CODE ------------------------------
