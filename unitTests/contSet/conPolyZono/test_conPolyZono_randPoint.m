function res = test_conPolyZono_randPoint
% test_conPolyZono_randPoint - unit test function of randPoint
%
% Syntax:
%    res = test_conPolyZono_randPoint
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

% Authors:       Tobias Ladner
% Written:       05-March-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% construct constrained polynomial zonotope
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
cPZ = conPolyZono(c,G,E,A,b,EC);

% check standard
p = randPoint(cPZ,10);
assert(all(contains(zonotope(cPZ),p)))

% check extreme
p = randPoint(cPZ,10,'extreme');
assert(all(contains(zonotope(cPZ),p)))

% check empty
p = randPoint(conPolyZono.empty(2));
assert(all(size(p) == [2 0]));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
