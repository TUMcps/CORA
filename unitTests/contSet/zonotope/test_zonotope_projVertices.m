function res = test_zonotope_projVertices
% test_zonotope_projVertices - unit test function for computation of
%    vertices of a 2D projection
%
% Syntax:
%    res = test_zonotope_projVertices
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

% 2D zonotope
c = [0;0];
G = [1.5 -1.5 0.5; 1 0.5 -1];
Z = conZonotope(c,G);

% compute vertices
V = vertices(Z);
V_proj = projVertices(Z);

% check vertices
if ~compareMatrices(V,V_proj,1e-14)
    res = false;
end


% 3D zonotope
c = [0; 0; 1];
G = [1.5 -1.5 0.5; 1 0.5 -1; 0 -0.5 -1];
Z = zonotope(c,G);

% compute vertices of full zonotope
V = vertices(Z);
% dimensions for projection
dims = {[1,2],[2,3],[1,3]};

% check all three projections
for i=1:length(dims)
    % computed vertices of projected zonotope
    V_proj = projVertices(Z,dims{i});

    % check vertices
    if ~compareMatrices(V_proj,V(dims{i},:),1e-14,'subset')
        res = false;
    end
end


% degenerate zonotope (line)
Z = zonotope([1;1],[1;-1]);

% compute vertices
V = [2 0; 0 2]';
V_proj = projVertices(Z);

% check vertices
if ~compareMatrices(V,V_proj,1e-14)
    res = false;
end

% degenerate zonotope (point)
Z = zonotope([1;1]);

% compute vertices
V = [1;1];
V_proj = projVertices(Z);

% check vertices
if ~compareMatrices(V,V_proj,1e-14)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
