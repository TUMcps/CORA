function res = testFlaky_polytope_vertices
% testFlaky_polytope_vertices - unit test function of vertices
%
% Syntax:
%    res = testFlaky_polytope_vertices
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
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% tolerance (similar to solver tolerance since we use linear programs)
tol = 1e-6;

% number of tests
nrTests = 50;

for i=1:nrTests

    % random dimension
    n = randi(5);

    % only n halfspaces -> unbounded (should throw an error)
    if n > 1
        % we can compute -+Inf vertices in 1D
        P = polytope(randn(n,n),ones(n,1));
        try
            V = vertices(P);
            throw(CORAerror('CORA:testFailed'));
        end
    end

    % sample semi-random vertices (no redundancies)
    V = -1+2*rand(n,2*n);
    V = V + [diag(10*ones(n,1)), -diag(10*ones(n,1))];

    % init polytope
    P = polytope(V);

    % compute vertices of instantiated polytope
    V_ = vertices(P);

    % check result
    if ~compareMatrices(V,V_,tol)
        throw(CORAerror('CORA:testFailed'));
    end


    % init random zonotope
    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',n+2);
    % delete aligned generators (to avoid numerical issues)
    Z = compact(Z,'aligned',1e-6);
    % lower bound for length of shortest generator with respect to longest
    % generator (to avoid numerical issues)
    ratio = vecnorm(Z.G) ./ max(vecnorm(Z.G));
    minratio = 0.01;
    if any(ratio < minratio)
        Z = zonotope(Z.c,max([minratio./ratio;ones(size(ratio))]).*Z.G);
    end

    % convert to polytope
    P = polytope(Z);

    % compute vertices
    V = vertices(Z);
    V_ = vertices(P);

    if ~compareMatrices(V,V_,tol)
        throw(CORAerror('CORA:testFailed'));
    end


    % degenerate polytopes: compute the vertices of a non-degenerate
    % polytope, add one dimension and compute the vertices again
    P = polytope.generateRandom('Dimension',n,'IsDegenerate',false);
    V = vertices(P);
    
    A = [P.A, zeros(length(P.b),1)]; b = P.b;
    Ae = [zeros(1,n) 1]; be = 0;
    P_ = polytope(A,b,Ae,be);
    V_ = vertices(P_);

    % compare vertices
    if ~compareMatrices(V,V_(1:n,:),tol)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
