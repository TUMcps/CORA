function res = test_polyZonotope_in
% test_polyZonotope_in - unit test function for containment checks of
%    polynomial zonotopes
%
% Syntax:  
%    res = test_polyZonotope_in
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

res = false;

% Analytical Test --------------------------------------------------------- 

% define polynomial zonotopes and points
pZ = polyZonotope([0;0],[2 0 1;0 2 1],[0.5;0],[1 0 3;0 1 1]);
p1 = [1;1];
p2 = [-1;3];
obj1 = polyZonotope([0;0],0.3*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);
obj2 = polyZonotope([0;0],0.4*[1 -2 1; 2 3 1],[],[1 0 2;0 1 1]);

% containment checks
res1 = in(pZ,p1,'approx');
res2 = in(pZ,p2,'approx');
res3 = in(pZ,obj1,'approx');
res4 = in(pZ,obj2,'approx');

% check if the result is correct
if res1 ~= 1 || res2 ~= 0 || res3 ~= 1 || res4 ~= 0
    error('test_polyZonotope_in: analytical test failed!');
end



res = true;

%------------- END OF CODE --------------