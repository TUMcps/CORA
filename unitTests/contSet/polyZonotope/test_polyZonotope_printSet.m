function res = test_polyZonotope_printSet
% test_polyZonotope_printSet - unit test function of printSet
%
% Syntax:
%    res = test_polyZonotope_printSet
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
pZ = polyZonotope.empty(2);

printSet(pZ)
printSet(pZ,'high')
printSet(pZ,'high',true)
printSet(pZ,'high',false)

% test normal set
c = [0;0];
G = [2 0 1;0 2 1];
GI = [0;0.5];
E = [1 0 3;0 1 1];

pZ = polyZonotope(c,G,GI,E);

printSet(pZ)
printSet(pZ,'high')
printSet(pZ,'high',true)
printSet(pZ,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
