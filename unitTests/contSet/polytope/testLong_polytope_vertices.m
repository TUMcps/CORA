function res = testLong_polytope_vertices
% testLong_polytope_vertices - unit test function of vertices
%
% Syntax:
%    res = testLong_polytope_vertices
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

% tolerance
tol = 1e-10;

% number of tests
nrTests = 50;

for i=1:nrTests

    % random dimension
    n = randi(5);

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
    
    Ae = [zeros(1,n) 1]; be = 0;
    P_ = polytope([P.A, zeros(length(P.b),1)],P.b,Ae,be);
    V_ = vertices(P_);

    % compare vertices
    if ~compareMatrices(V,V_(1:n,:),tol)
        throw(CORAerror('CORA:testFailed'));
    end

end

% ------------------------------ END OF CODE ------------------------------
