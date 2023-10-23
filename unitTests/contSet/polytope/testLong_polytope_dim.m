function res = testLong_polytope_dim
% testLong_polytope_dim - unit test function for dim
%
% Syntax:
%    res = testLong_polytope_dim()
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests
    
    % random dimension
    n = randi(5);

    % instantiation: generateRandom
    P = polytope.generateRandom('Dimension',n);

    % check dimension
    if dim(P) ~= n
        res = false; return
    end

    % instantiation: A and b matrix
    nrCon = randi(2*n);
    A = randn(nrCon,n);
    b = randn(nrCon,1);

    P = polytope(A,b);

    % check dimension
    if dim(P) ~= n
        res = false; return
    end

    % instantiation: vertices
    V = randn(n,randi([n+1,2*n]));

    P = polytope(V);

    % check dimension
    if dim(P) ~= n
        res = false; return
    end

    % instantiation: equality representation
    nrCon = randi(2*n);
    Ae = randn(nrCon,n);
    be = randn(nrCon,1);

    P = polytope([],[],Ae,be);

    % check dimension
    if dim(P) ~= n
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
