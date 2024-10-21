function res = test_conPolyZono_printSet
% test_conPolyZono_printSet - unit test function of printSet
%
% Syntax:
%    res = test_conPolyZono_printSet
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

% Authors:       Tobias Ladner
% Written:       10-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% test empty
cPZ = conPolyZono.empty(2);

printSet(cPZ)
printSet(cPZ,'high')
printSet(cPZ,'high',true)
printSet(cPZ,'high',false)

% test normal set
c = [0;0];
G = [1 0 1 -1; 0 1 1 1];
E = [1 0 1 2; 0 1 1 0; 0 0 1 1];
A = [1 -0.5 0.5];
b = 0.5;
EC = [0 1 2; 1 0 0; 0 1 0];
cPZ = conPolyZono(c,G,E,A,b,EC);

printSet(cPZ)
printSet(cPZ,'high')
printSet(cPZ,'high',true)
printSet(cPZ,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
