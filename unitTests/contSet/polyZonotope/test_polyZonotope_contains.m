function res = test_polyZonotope_contains
% test_polyZonotope_contains - unit test function for containment checks of
%    polynomial zonotopes
%
% Syntax:
%    res = test_polyZonotope_contains
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

% Authors:       Niklas Kochdumper
% Written:       13-January-2020
% Last update:   10-December-2023 (added test for case that prev. failed)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% assume true
res = true(0);

% define polynomial zonotopes and points
pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0.5;0],[1 0 3;0 1 1]);
p1 = [1;1];
p2 = [-1;3];
pZ1 = polyZonotope([0;0],0.3*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
pZ2 = polyZonotope([0;0],0.4*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
box = polyZonotope([0;0],[1 0; 0.0 1],[],[1 0;0 1]);
smallBox = 0.5*box;

% containment checks
res(end+1,1) = contains(pZ,p1,'approx');
res(end+1,1) = ~contains(pZ,p2,'approx');
res(end+1,1) = contains(pZ,pZ1,'approx');
res(end+1,1) = ~contains(pZ,pZ2,'approx');
res(end+1,1) = contains(box,smallBox,'approx');

% empty set
% pZ_empty = polyZonotope.empty(2);
% res(end+1,1) = ~contains(pZ_empty,pZ1) && contains(pZ1,pZ_empty);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
