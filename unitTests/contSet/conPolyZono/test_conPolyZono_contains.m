function res = test_conPolyZono_contains
% test_conPolyZono_contains - unit test function of contains
%
% Syntax:
%    res = test_conPolyZono_contains
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
% See also: testLong_conPolyZono_contains

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

I = interval([-1;-1],[1;1])*0.1;

% check containment 'exact'
assertThrowsAs(@contains,'CORA:noExactAlg',cPZ,c);
assertThrowsAs(@contains,'CORA:noExactAlg',cPZ,I);

% check containment 'approx' (just if it runs through)
contains(cPZ,c,'approx');

% other checks in long test ---

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
