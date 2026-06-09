function res = testLong_polytope_projVertices
% testLong_polytope_projVertices - unit test function for
%    computation of vertices of a 2D projection
%
% Syntax:
%    res = testLong_polytope_projVertices
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

% Authors:       Niklas Kochdumper
% Written:       20-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i = 1:nrTests
    
    % random dimension
    n = randi([2,4]);

    % random polytope
    P = polytope.generateRandom("Dimension",n);
    
    % random dimensions for projection
    projDims = randperm(n,2);

    % compute vertices
    V_proj = projVertices(P,projDims);
    V = vertices(P);
    
    % check vertices
    assertLoop(compareMatrices(V_proj,V(projDims,:),1e-10,'subset'),i)
end

% ------------------------------ END OF CODE ------------------------------
