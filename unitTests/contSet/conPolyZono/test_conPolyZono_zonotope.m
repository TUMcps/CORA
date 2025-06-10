function res = test_conPolyZono_zonotope
% test_conPolyZono_zonotope - unit test function of constructor
%
% Syntax:
%    res = test_conPolyZono_zonotope
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
% Written:       22-May-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% very simple test, majority is done in the long test

% init cPZ
c = [0;0];
G = [2 1 2 1; 0 2 2 1];
E = [1 0 2 0; 0 1 1 0; 0 0 0 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [1 0 0; 0 1 2; 0 1 0];
cPZ = conPolyZono(c,G,E,A,b,EC);

% init true Z
c = [0; 0];
G = [ 2 1 2 1 ; 0 2 2 1 ];
Z_true = zonotope(c,G);

% compute zonotope enclosure
Z = zonotope(cPZ);

% check equality
assert(isequal(Z,Z_true))

% test completed
res = true;

end

% ------------------------------ END OF CODE ------------------------------
