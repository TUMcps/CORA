function res = test_capsule_printSet
% test_capsule_printSet - unit test function of printSet
%
% Syntax:
%    res = test_capsule_printSet
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
C = capsule.empty(2);

printSet(C)
printSet(C,'high')
printSet(C,'high',true)
printSet(C,'high',false)

% test normal set
c = [1;2];
g = [2;1];
r = 1;
C = capsule(c,g,r);

printSet(C)
printSet(C,'high')
printSet(C,'high',true)
printSet(C,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
