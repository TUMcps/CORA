function res = test_ellipsoid_vertices
% test_ellipsoid_vertices - unit test function of vertices
%
% Syntax:
%    res = test_ellipsoid_vertices
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
% Written:       25-July-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
n = 2;
E = ellipsoid.empty(n);
V = vertices(E);
res(end+1,1) = isempty(V) && isnumeric(V) && size(V,1) == n;

% 1D case, just a point
E = ellipsoid(0,1);
V = vertices(E);
res(end+1,1) = size(V,2) == 1 && withinTol(V,1);

% 1D case, bounded line
E = ellipsoid(4,-1);
V = vertices(E);
res(end+1,1) = size(V,2) == 2 && compareMatrices(V,[-3,1]);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
