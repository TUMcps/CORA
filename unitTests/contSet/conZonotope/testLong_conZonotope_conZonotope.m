function res = testLong_conZonotope_conZonotope
% testLong_conZonotope_conZonotope - unit test function of
%    conZonotope (constructor)
%
% Syntax:
%    res = testLong_conZonotope_conZonotope
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
% Written:       19-March-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;
nrOfTests = 1000;
for i=1:nrOfTests

    % random dimension
    n = randi(25);
    % number of generators
    nrGens = randi(25);
    % number of constraints
    nrCon = randi(10);
    
    % random center, generator matrix (-> full matrix of zonotope),
    % as well as constraint matrix, constraint vector
    c = randn(n,1);
    G = randn(n,nrGens);
    Z = [c,G];
    A = randn(nrCon,nrGens);
    b = randn(nrCon,1);
    
    % admissible initializations
    % only zonotope matrix
    cZ = conZonotope(Z);
    assertLoop(compareMatrices([cZ.c,cZ.G],Z),i)

    % center and generator matrix
    cZ = conZonotope(c,G);
    assertLoop(all(all(withinTol([cZ.c,cZ.G],Z))),i)

    % zonotope matrix, constraint matrix, constraint vector
    cZ = conZonotope(Z,A,b);
    assertLoop(all(all(withinTol([cZ.c,cZ.G],Z))),i)
    assertLoop(all(all(withinTol(cZ.A,A))),i)
    assertLoop(all(withinTol(cZ.b,b)),i)
    
    % center, generator matrix, constraint matrix, constraint vector
    cZ = conZonotope(c,G,A,b);
    assertLoop(all(all(withinTol([cZ.c,cZ.G],Z))),i)
    assertLoop(all(all(withinTol(cZ.A,A))),i)
    assertLoop(all(withinTol(cZ.b,b)),i)
    
    
    % wrong initializations
    c_plus1 = randn(n+1,1);
    c_mat = randn(n);
    G_plus1 = randn(n,nrGens+1);
    Z_plus1 = [c,G_plus1];
    A_plus1 = randn(nrCon,nrGens+1);
    b_plus1 = randn(nrCon+1,1);
    
    % center and generator matrix of different dimensions
    assertThrowsAs(@conZonotope,'CORA:wrongInputInConstructor',c_plus1,G);
    
    % A does not fit Z
    assertThrowsAs(@conZonotope,'CORA:wrongInputInConstructor',Z,A_plus1,b);
    assertThrowsAs(@conZonotope,'CORA:wrongInputInConstructor',Z_plus1,A,b);
    
    % A does not fit b
    assertThrowsAs(@conZonotope,'CORA:wrongInputInConstructor',Z,A,b_plus1);
    
    % center is a matrix
    if n ~= 1
        assertThrowsAs(@conZonotope,'CORA:wrongInputInConstructor',c_mat,G);
    end
    
    % too many input arguments
    assertThrowsAs(@conZonotope,'CORA:numInputArgsConstructor',c,G,A,b,b);
end

% ------------------------------ END OF CODE ------------------------------
