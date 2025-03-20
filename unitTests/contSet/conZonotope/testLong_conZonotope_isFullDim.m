function res = testLong_conZonotope_isFullDim
% testLong_conZonotope_isFullDim - unit test function of isFullDim
%
% Syntax:
%    res = testLong_conZonotope_isFullDim
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
% Written:       14-March-2021
% Last update:   04-February-2024 (AK, added degenerate cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% number of tests
nrOfTests = 1000;

for i=1:nrOfTests
    % random dimension
    n = randi([1,50]);
    % random center
    c = randn(n,1);
    % random generator matrix (at least n generators)
    G = randn(n,n+randi(10));
    nrGens = size(G,2);
    
    % instantiate conZonotope without constraints
    cZ = conZonotope(c,G);
    
    % assert correctness
    assertLoop(isFullDim(cZ),i)
    
    % random constraints so that conZonotope represents just a point
    % as A being diagonal forces each independent factor to one value
    A = diag(1+abs(rand(nrGens,1)));
    b = sign(randn(nrGens,1));
    % instantiate conZonotope with constraints
    cZ = conZonotope(c,G,A,b);
    
    % assert correctness
    [res_subspace, subspace] = isFullDim(cZ);
    assert(~isFullDim(cZ));
    assert(~res_subspace);
    assert(isempty(subspace));
    
end

for i=1:(nrOfTests/10)
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

    % Now dealing with the constraints:
    % random number of constraints
    n_constraints = randi([1,n+1]);

    b = randn([n_constraints 1]);
    A = randn([n_constraints gamma]);
    
    % instantiate zonotope
    cZ = conZonotope(c,Gdeg, A, b);
    
    % check if full-dimensional
    res_single = isFullDim(cZ);
    [res_subspace, subspace] = isFullDim(cZ);

    % transform into polytope, to check if both methods yield the same
    P = polytope(cZ);
    [res_poly,P_subspace] = isFullDim(P);

    same_subspace = rank([subspace P_subspace], 1e-5) == rank(subspace);

    assertLoop(res_single == res_poly,i);
    assertLoop(res_subspace == res_poly,i);
    assertLoop(same_subspace,i);
end

% ------------------------------ END OF CODE ------------------------------
