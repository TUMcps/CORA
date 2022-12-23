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
%    res - boolean 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Niklas Kochdumper
% Written:      13-January-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% assume true
res = true;

% define polynomial zonotopes and points
pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0.5;0],[1 0 3;0 1 1]);
p1 = [1;1];
p2 = [-1;3];
obj1 = polyZonotope([0;0],0.3*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
obj2 = polyZonotope([0;0],0.4*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);

% empty set
obj_e = polyZonotope();
res_e = ~contains(obj_e,obj1) && contains(obj1,obj_e);

% containment checks
res1 = contains(pZ,p1,'approx');
res2 = contains(pZ,p2,'approx');
res3 = contains(pZ,obj1,'approx');
res4 = contains(pZ,obj2,'approx');

% check if the result is correct
if ~res1 || res2 || ~res3 || res4 || ~res_e
    res = false;
end

%------------- END OF CODE --------------