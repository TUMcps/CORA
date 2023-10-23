function res = testLong_interval_projVertices
% testLong_interval_projVertices - unit test function for
%    computation of vertices of a 2D projection
%
% Syntax:
%    res = testLong_interval_projVertices
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

% Authors:       Mark Wetzlinger
% Written:       21-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true;

% number of tests
nrTests = 100;

for i=1:nrTests
    
    % random dimension
    n = randi([2,10]);

    % instantiate randomized zonotope
    I = interval.generateRandom('Dimension',n);

    % random projection
    projDims = randperm(n,2);

    % compute vertices of interval
    V = vertices(I);

    % compute vertices of projected interval
    V_ = vertices(project(I,projDims));

    % compute projected vertices
    V_proj = projVertices(I,projDims);

    % has to have 4 vertices (all dimensions have non-zero extension)
    if size(V_proj,2) ~= 4
        res = false; return
    end

    % check vertices
    if ~compareMatrices(V_proj,V(projDims,:),1e-14,'subset') ...
            || ~compareMatrices(V_,V_proj,1e-14)
        res = false; return
    end

end

% ------------------------------ END OF CODE ------------------------------
