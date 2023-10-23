function res = testLong_polytope_conZonotope
% testLong_polytope_conZonotope - unit test function for conversion
%    to conZonotope objects
%
% Syntax:
%    res = testLong_polytope_conZonotope
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
% Written:       04-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 10;
% tolerance
tol = 1e-8;

for i=1:nrTests

    % random dimension
    n = randi(6);

    % instantiation: generateRandom
    P = polytope.generateRandom('Dimension',n,'NrConstraints',2*n);
    % compute vertices
    V_P = vertices(P);

    % convert to constrained zonotope
    cZ_sF = conZonotope(P,'exact:supportFunc');
    cZ_vert = conZonotope(P,'exact:vertices');
    

    % compute vertices
    V_sF = vertices(cZ_sF);
    V_vert = vertices(cZ_vert);

    % compare vertices
    if ~compareMatrices(V_P,V_sF,tol) || ~compareMatrices(V_P,V_vert,tol)
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
