function res = test_conZonotope_projVertices
% test_conZonotope_projVertices - unit test function for computation of
%    vertices of a 2D projection
%
% Syntax:
%    res = test_conZonotope_projVertices
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

% 2D constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1]; b = 1;
cZ = conZonotope(Z,A,b);

% compute vertices
V = vertices(cZ);
V_proj = projVertices(cZ);

% check vertices
if ~compareMatrices(V,V_proj,1e-14)
    res = false;
end


% higher-dimensional constrained zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1; 1 0 -0.5 -1];
A = [1 1 1]; b = 1;
cZ = conZonotope(Z,A,b);

% compute vertices of full constrained zonotope
V = vertices(cZ);
% dimensions for projection
dims = {[1,2],[2,3],[1,3]};

% check all three projections
for i=1:length(dims)
    % computed vertices of projected constrained zonotope
    V_proj = projVertices(cZ,dims{i});

    % check vertices
    if ~compareMatrices(V_proj,V(dims{i},:),1e-14,'subset')
        res = false;
    end
end


% convex hull of two constrained zonotopes
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1]; b = 1;
cZ1 = conZonotope(Z,A,b);

Z = [4 2 0 0;4 1 1 0];
A = [1 1 -1]; b = 0;
cZ2 = conZonotope(Z,A,b);

% compute convex hull
cZ_ = convHull(cZ1,cZ2);

% compute vertices
V = vertices(cZ_);
V_proj = projVertices(cZ_);

% check vertices
if ~compareMatrices(V_proj,V,1e-14,'subset')
    res = false;
end


% ------------------------------ END OF CODE ------------------------------
