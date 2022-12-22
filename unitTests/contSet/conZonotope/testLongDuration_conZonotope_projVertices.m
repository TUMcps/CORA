function res = testLongDuration_conZonotope_projVertices
% testLongDuration_conZonotope_projVertices - unit test function for
%    computation of vertices of a 2D projection
%
% Syntax:  
%    res = testLongDuration_conZonotope_projVertices
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

    % random center
    c = randn(n,1);
    % random generator matrix
    G = randn(n,n+floor(rand*n));

    % constraint such that factor space not empty
    A = randn(1,size(G,2));
    b = 0;

    % instantiate constrained zonotope
    cZ = conZonotope(c,G,A,b);
    
    % random dimensions for projection
    projDims = randperm(n,2);

    % compute vertices
    V = vertices(project(cZ,projDims));
    V_proj = projVertices(cZ,projDims);
    
    % check vertices
    if ~compareMatrices(V,V_proj,1e-14)
        res = false;
    end

end

%------------- END OF CODE --------------