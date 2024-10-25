function res = test_capsule_vertices
% test_capsule_vertices - unit test function of polygon
%
% Syntax:
%    res = test_capsule_vertices
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       24-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init 2D capsule
C = capsule([1;-1],[4;2],1);

% compute polygon (inner-approximation)
V = vertices(C);

% check whether all points are contained in the capsule
assert(all(contains(C,V)));

% check only center
c = [1;2];
C = capsule(c);
assert(compareMatrices(c,vertices(C)));

res = true;

% ------------------------------ END OF CODE ------------------------------
