function res = test_polygon_isequal()
% test_polygon_isequal - unit test function for equality check of polygons
%
% Syntax:
%    res = test_polygon_isequal
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
% See also: polygon

% Authors:       Mark Wetzlinger
% Written:       27-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = [];

% init polygon via vertices
V = [4 0; 3 2; 0 3; -3 0; -2 -3; 0 -3; 2 -2]';
P1 = polygon(V);
% has to be equal to itself
res(end+1,1) = isequal(P1,P1);

% re-order points
V = [V(:,3:7) V(:,1:2)];
P2 = polygon(V);
% still the same shape
res(end+1,1) = isequal(P1,P2);

% different polygon by a little
V_ = [4 0; 3 2; 0 3; -3 1e-6; -2 -3; -1e-6 -3; 2 -2]';
P2 = polygon(V_);
% not the same using standard tolerance
res(end+1,1) = ~isequal(P1,P2);
% indeed the same using large enough tolerance
res(end+1,1) = isequal(P1,P2,1e-4);

% polygon with different number of vertices
V_ = [4 0; 3 2; 0 3; -3 0; 0 -3; 2 -2]';
P2 = polygon(V_);
% not the same shape
res(end+1,1) = ~isequal(P1,P2);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
