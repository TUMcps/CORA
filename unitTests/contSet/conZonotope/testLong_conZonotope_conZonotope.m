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
    if ~compareMatrices([cZ.c,cZ.G],Z)
        res = false; break;
    end

    % center and generator matrix
    cZ = conZonotope(c,G);
    if ~all(all(withinTol([cZ.c,cZ.G],Z)))
        res = false; break;
    end

    % zonotope matrix, constraint matrix, constraint vector
    cZ = conZonotope(Z,A,b);
    if ~all(all(withinTol([cZ.c,cZ.G],Z))) || ~all(all(withinTol(cZ.A,A))) ...
            || ~all(withinTol(cZ.b,b))
        res = false; break;
    end
    
    % center, generator matrix, constraint matrix, constraint vector
    cZ = conZonotope(c,G,A,b);
    if ~all(all(withinTol([cZ.c,cZ.G],Z))) || ~all(all(withinTol(cZ.A,A))) ...
            || ~all(withinTol(cZ.b,b))
        res = false; break;
    end
    
    
    % wrong initializations
    c_plus1 = randn(n+1,1);
    c_mat = randn(n);
    G_plus1 = randn(n,nrGens+1);
    Z_plus1 = [c,G_plus1];
    A_plus1 = randn(nrCon,nrGens+1);
    b_plus1 = randn(nrCon+1,1);
    
    % center and generator matrix of different dimensions
    try
        cZ = conZonotope(c_plus1,G); % <- should throw error here
        res = false; break;
    end
    
    % A does not fit Z
    try
        cZ = conZonotope(Z,A_plus1,b); % <- should throw error here
        res = false; break;
    end
    try
        cZ = conZonotope(Z_plus1,A,b); % <- should throw error here
        res = false; break;
    end
    
    % A does not fit b
    try
        cZ = conZonotope(Z,A,b_plus1); % <- should throw error here
        res = false; break;
    end
    
    % center is a matrix
    if n ~= 1
        try
            cZ = conZonotope(c_mat,G); % <- should throw error here
            res = false; break;
        end
    end
    
    % too many input arguments
    try
        cZ = conZonotope(c,G,A,b,b); % <- should throw error here
        res = false; break;
    end
end

% ------------------------------ END OF CODE ------------------------------
