function res = test_polytope_printSet
% test_polytope_printSet - unit test function of printSet
%
% Syntax:
%    res = test_polytope_printSet
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
P = polytope.empty(2);

printSet(P)
printSet(P,'high')
printSet(P,'high',true)
printSet(P,'high',false)

% test normal set in Hrep
A = [1 0 -1 0 1; 0 1 0 -1 1]';
b = [3; 2; 3; 2; 1];
P = polytope(A,b);

printSet(P)
printSet(P,'high')
printSet(P,'high',true)
printSet(P,'high',false)

% test normal set in Vrep
V = [1 0 1; 0 1 1];
P = polytope(V);

printSet(P)
printSet(P,'high')
printSet(P,'high',true)
printSet(P,'high',false)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
