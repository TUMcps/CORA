function res = test_contSet_projVertices
% test_contSet_projVertices - regression test for the index-shift bug in
%    aux_projVerticesSupportFunc (see @contSet/projVertices.m).
%
% Syntax:
%    res = test_contSet_projVertices
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
% See also: projVertices

% Authors:       Tobias Ladner
% Written:       19-February-2026
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% found edge cases using Claude Code:

% --- minimal index-shift regression case ---------------------------------
%
% Construct a 3D conZonotope whose 2D projection to dims [1,2] is a
% quadrilateral.  The generators and constraint are chosen so that the
% support at the 120-degree direction hits the midpoint of a projected edge
% rather than a vertex, making that initial vertex collinear with its
% neighbour (the vertex inserted by the main refinement loop).  In the
% buggy code this leads to idx3 pointing beyond the end of V after the
% collinear vertex is removed.

% generators (columns) embedded in 3D, center at origin
G = [1  0  1  0; ...
     0  1  0  1; ...
     1 -1  0  0];
c = [0; 0; 0];

% equality constraint: forces the set to have a face perpendicular to the
% 120-degree support direction when projected to the [1,2] plane
A_eq = [1, 1, -1, -1];
b_eq = 0;

cZ = conZonotope([c, G], A_eq, b_eq);

% supportFunc algorithm: would crash (index out of bounds) with old code
% if the LP returns the midpoint of a projected edge as an initial vertex
V_sf = projVertices(cZ, [1,2], 'supportFunc');

% verify the result is a subset of the angle-algorithm result (or equal)
V_ang = projVertices(cZ, [1,2], 'angle');
assert(compareMatrices(V_sf, V_ang, 1e-6, 'subset') || ...
       compareMatrices(V_ang, V_sf, 1e-6, 'subset'), ...
       "supportFunc and angle algorithms disagree on projection vertices.");


% --- second regression case: multiple initial collinear vertices ---------
%
% A conZonotope in 4D projected to dims [1,2].  The equality constraint
% links generators so that both the 120-degree and 240-degree initial
% support vectors can land on an edge interior, making both idx2 and idx3
% candidates for deletion.  The old code would either go out of bounds or
% silently delete the wrong vertex.

G2 = [2  0  0  1  1; ...
      0  2  0  1 -1; ...
      0  0  2 -1  0; ...
      0  0  0  0  1];
c2 = zeros(4, 1);
A_eq2 = [1 -1  0  0  0; ...
          0  0  1 -1  0];
b_eq2 = [0; 0];

cZ2 = conZonotope([c2, G2], A_eq2, b_eq2);

V_sf2  = projVertices(cZ2, [1,2], 'supportFunc');
V_ang2 = projVertices(cZ2, [1,2], 'angle');
assert(compareMatrices(V_sf2, V_ang2, 1e-6, 'subset') || ...
       compareMatrices(V_ang2, V_sf2, 1e-6, 'subset'), ...
       "supportFunc and angle algorithms disagree (second case).");


% --- basic sanity check: 2D zonotope (no deletion expected) --------------

Z = zonotope([0;0], [1 0; 0 1]);
V_sf3  = projVertices(Z, [1,2], 'supportFunc');
V_ang3 = projVertices(Z, [1,2], 'angle');
assert(compareMatrices(V_sf3, V_ang3, 1e-6, 'subset') || ...
       compareMatrices(V_ang3, V_sf3, 1e-6, 'subset'), ...
       "supportFunc and angle algorithms disagree for zonotope.");


% --- idx2/idx3 empty: no deletions at duplicate-check stage --------------
%
% When the 120-degree or 240-degree support returns a duplicate of v1, the
% duplicate is removed early and idx2 or idx3 is set to [].  The batch
% code must handle empty indices (they are dropped by MATLAB concatenation
% [1, [], idx3] -> [1, idx3]).

Z2 = zonotope([0;0], [1 0; 0 0]);   % degenerate: line segment along x
V_sf4 = projVertices(Z2, [1,2], 'supportFunc');
assert(size(V_sf4, 1) == 2, "Output must have 2 rows.");
assert(size(V_sf4, 2) >= 1, "Output must have at least 1 vertex.");

end

% ------------------------------ END OF CODE ------------------------------
