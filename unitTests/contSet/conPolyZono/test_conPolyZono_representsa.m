function res = test_conPolyZono_representsa
% test_conPolyZono_representsa - unit test function of representation check
%
% Syntax:
%    res = test_conPolyZono_representsa
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

% Authors:       Niklas Kochdumper
% Written:       04-February-2021
% Last update:   10-January-2024 (MW, copied from removed isempty function)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% instantiate constrained polynomial zonotope
c = [0;0];
G = [1 0 1;0 1 1];
E = [1 0 2;0 1 1];
A = [1 -1 0; 0 -1 1];
b1 = [0; 1];
EC = [2 0 1; 0 1 0];
cPZ1 = conPolyZono(c,G,E,A,b1,EC);
% b2 = [0; 0];
% cPZ2 = conPolyZono(c,G,E,A,b2,EC);

res(end+1,1) = representsa(cPZ1,'emptySet',1e-8,'linearize',3,7);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
