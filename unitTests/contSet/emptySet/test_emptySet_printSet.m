function res = test_emptySet_printSet
% test_emptySet_printSet - unit test function of printSet
%
% Syntax:
%    res = test_emptySet_printSet
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
O = emptySet.empty(2);

printSet(O)
printSet(O,'high')
printSet(O,'high',true)
printSet(O,'high',false)

% test normal set
n = 2;
O = emptySet(n);

printSet(O)
printSet(O,'high')
printSet(O,'high',true)
printSet(O,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
