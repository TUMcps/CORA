function res = test_ellipsoid_printSet
% test_ellipsoid_printSet - unit test function of printSet
%
% Syntax:
%    res = test_ellipsoid_printSet
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
E = ellipsoid.empty(2);

printSet(E)
printSet(E,'high')
printSet(E,'high',true)
printSet(E,'high',false)

% test normal set
Q = [2.7 -0.2;-0.2 2.4];
q = [1;2];
E = ellipsoid(Q, q);

printSet(E)
printSet(E,'high')
printSet(E,'high',true)
printSet(E,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
