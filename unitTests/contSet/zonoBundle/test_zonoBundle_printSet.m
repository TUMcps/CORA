function res = test_zonoBundle_printSet
% test_zonoBundle_printSet - unit test function of printSet
%
% Syntax:
%    res = test_zonoBundle_printSet
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
zB = zonoBundle.empty(2);

printSet(zB)
printSet(zB,'high')
printSet(zB,'high',true)
printSet(zB,'high',false)

% test normal set
Z1 = zonotope([1 3 0; 1 0 2]);
Z2 = zonotope([0 2 2; 0 2 -2]);
zB = zonoBundle({Z1,Z2});

printSet(zB)
printSet(zB,'high')
printSet(zB,'high',true)
printSet(zB,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
