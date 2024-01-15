function res = test_conPolyZono_conPolyZono
% test_conPolyZono_conPolyZono - unit test function of constructor
%
% Syntax:
%    res = test_conPolyZono_conPolyZono
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
% Written:       15-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test example in doctring
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];

cPZ = conPolyZono(c,G,E,A,b,EC);

% test variants in syntax of docstring
GI = [4 1; 0 2];
id = [1;2;4];

cPZ = conPolyZono(c,G,E);
cPZ = conPolyZono(c,G,E,GI);
cPZ = conPolyZono(c,G,E,GI,id);
cPZ = conPolyZono(c,G,E,A,b,EC);
cPZ = conPolyZono(c,G,E,A,b,EC,GI);
cPZ = conPolyZono(c,G,E,A,b,EC,GI,id);

% test copy constructor
cPZ = conPolyZono(cPZ);

% test empty set
cPZ = conPolyZono.empty(2);

res = true;

% ------------------------------ END OF CODE ------------------------------
