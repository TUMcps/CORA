function res = test_ellipsoid_generators
% test_ellipsoid_generators - unit test function of generators
%
% Syntax:
%    res = test_ellipsoid_generators
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

% Authors:       Adrian Kulmburg
% Written:       05-March-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty ellipsoid
E = ellipsoid.empty(2);
G = generators(E);
res = isempty(G) && isnumeric(G) && all(size(G) == [2,0]);

% 2D zonotope
c = [-2; 1];
G = [1 2 0;2 3 1];
E = G * ellipsoid(eye(3)) + c;
G_ = generators(E);
E_ = G_ * ellipsoid(eye(size(G_, 2))) + c;
res(end+1,1) = E == E_;

% degenerate ellipsoid
c = [-2; 1; 3];
G = [1 2 0 5;2 3 1 -1;0 0 0 0];
E = G * ellipsoid(eye(4)) + c;
G_ = generators(E);
E_ = G_ * ellipsoid(eye(size(G_, 2))) + c;
res(end+1,1) = E == E_;

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
