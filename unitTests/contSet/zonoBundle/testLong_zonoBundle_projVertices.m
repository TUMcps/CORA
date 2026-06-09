function res = testLong_zonoBundle_projVertices
% testLong_zonoBundle_projVertices - unit test function for
%    computation of vertices of a 2D projection
%
% Syntax:
%    res = testLong_zonoBundle_projVertices
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
nrTests = 5;

for i = 1:nrTests
    
    % random dimension and number of zonotopes
    n = randi([2,4]);
    nrZonos = randi([2,4]);

    % random zonotope bundle
    zB = zonoBundle.generateRandom('Dimension',n,'NrZonotopes',nrZonos);
    
    % random dimensions for projection
    projDims = randperm(n,2);

    % compute projected vertices
    V_proj1 = projVertices(zB,projDims,'angle');
    V_proj2 = projVertices(zB,projDims,'supportFunc');
    
    % check vertices (large tolerance since linear programs are used)
    assertLoop(compareMatrices(V_proj2,V_proj1,1e-5,'subset'),i)
end

% ------------------------------ END OF CODE ------------------------------
