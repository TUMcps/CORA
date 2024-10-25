function res = test_polytope_project
% test_polytope_project - unit test function for projection of polytopes
%
% Syntax:
%    res = test_polytope_project()
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

% Authors:       Niklas Kochdumper, Mark Wetzlinger
% Written:       21-December-2020
% Last update:   14-July-2024 (MW, test vertex representation)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 1D, projection to itself
A = [1; -1]; b = [1; 0];
P = polytope(A,b);
P_proj = project(P,1);
assert(P == P_proj);


% 2D, bounded
A = [-1 0; 2 4; 1 -2];
b = [-1; 14; -1];
P = polytope(A,b);
% project to dimension 1 and compare to true result
P_proj = project(P,1);
V = vertices(P_proj);
V_ = [1 3];
assert(compareMatrices(V,V_));
% project to dimension 2 and compare to true result
P_proj = project(P,2);
V = vertices(P_proj);
V_ = [1 3];
assert(compareMatrices(V,V_));

% 2D, unbounded
P = polytope([1 0; -1 0; 0 1], [1;1;1]);
% project to dimension 1
P_proj = project(P,1);
% check if it is still unbounded
assert(isBounded(P_proj));

% 2D, degenerate
P = polytope([1 0; -1 0; 0 1; 0 -1], [1;1;1;-1]);
P_proj = project(P,1);
P_true = polytope([1;-1], [1;1]);
% polytope should remain unchanged but now full dimensional dimension
assert(P_true == P_proj && isFullDim(P_proj));

% 2D, scaled unit square, vertex instantiation
V = [2 0; -2 0; 0 2; 0 -2]';
P = polytope(V);
% project to dimension 2 and compare to true result
P_proj = project(P,2);
V = vertices(P_proj);
V_ = [2 -2];
assert(compareMatrices(V,V_));


% 3D, fully empty
A = zeros(0,3); b = zeros(0,0);
P = polytope(A,b);
P_proj = project(P,1);
A_true = zeros(0,1); b_true = zeros(0,0);
P_true = polytope(A_true,b_true);
assert(isequal(P_proj,P_true));

% 3D, vertex representation
V = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1]';
P = polytope(V);
P_proj = project(P,[1,3]);
V_proj = vertices(P_proj);
V_true = [-1 1; 1 1; 1 -1]';
assert(compareMatrices(V_proj,V_true,1e-8,"equal",false));


% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
