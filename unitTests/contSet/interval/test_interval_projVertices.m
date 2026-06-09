function res = test_interval_projVertices
% test_interval_projVertices - unit test function for computation of
%    vertices of a 2D projection
%
% Syntax:
%    res = test_interval_projVertices
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

% instantiate interval
I = interval([-3;-4;-2],[1;4;3]);

% dimensions for projection
projDims = [2,3];

% true vertices
V = [-4 -2; -4 3; 4 3; 4 -2]';

% compute projected vertices
V_proj1 = projVertices(I,projDims,'angle');
V_proj2 = projVertices(I,projDims,'supportFunc');

% check vertices
assert(compareMatrices(V,V_proj1,1e-14))
assert(compareMatrices(V,V_proj2,1e-14))


% degenerate interval (line)
I = interval([-3;-4;3],[1;4;3]);

% dimensions for projection
projDims = [1,3];

% true vertices
V = [-3 3; 1 3]';

% compute projected vertices
V_proj1 = projVertices(I,projDims,'angle');
V_proj2 = projVertices(I,projDims,'supportFunc');

% check vertices
assert(compareMatrices(V,V_proj1,1e-14))
assert(compareMatrices(V,V_proj2,1e-14))


% degenerate interval (point)
I = interval([-3;-4;3],[1;-4;3]);

% dimensions for projection
projDims = [2,3];

% true vertices
V = [-4;3];

% compute projected vertices
V_proj1 = projVertices(I,projDims,'angle');
V_proj2 = projVertices(I,projDims,'supportFunc');

% check vertices
assert(compareMatrices(V,V_proj1,1e-14))
assert(compareMatrices(V,V_proj2,1e-14))

end

% ------------------------------ END OF CODE ------------------------------
