function res = testLongDuration_zonotope_projVertices
% testLongDuration_zonotope_projVertices - unit test function for
%    computation of vertices of a 2D projection
%
% Syntax:  
%    res = testLongDuration_zonotope_projVertices
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
% See also: none

% Author:       Mark Wetzlinger
% Written:      21-December-2022
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% number of tests
nrTests = 25;

for i=1:nrTests
    
    % random dimension
    n = randi([2,10]);

    % instantiate randomized zonotope
    Z = zonotope.generateRandom('Dimension',n,'NrGenerators',n+floor(n*rand));

    % random projection
    projDims = randperm(n,2);

    % compute vertices
    V = vertices(project(Z,projDims));

    % compute vertices in projection
    V_proj = projVertices(Z,projDims);

    % visualization
%     figure; hold on; box on;
%     plot(project(Z,projDims));
%     scatter(V(1,:),V(2,:),16,'r','filled');
%     scatter(V_proj(1,:),V_proj(2,:),16,'g');

    % check vertices
    if ~compareMatrices(V_proj,V,1e-14)
        res = false; return
    end

end

%------------- END OF CODE --------------
