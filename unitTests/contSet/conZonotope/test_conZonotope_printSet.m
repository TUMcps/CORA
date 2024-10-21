function res = test_conZonotope_printSet
% test_conZonotope_printSet - unit test function of printSet
%
% Syntax:
%    res = test_conZonotope_printSet
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
cZ = conZonotope.empty(2);

printSet(cZ)
printSet(cZ,'high')
printSet(cZ,'high',true)
printSet(cZ,'high',false)

% test normal set
Z = [0 3 0 1;0 0 2 1];
A = [1 0 1]; b = 1;
cZ = conZonotope(Z,A,b);

printSet(cZ)
printSet(cZ,'high')
printSet(cZ,'high',true)
printSet(cZ,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
