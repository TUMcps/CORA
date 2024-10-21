function res = test_interval_printSet
% test_interval_printSet - unit test function of printSet
%
% Syntax:
%    res = test_interval_printSet
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
I = interval.empty(2);

printSet(I)
printSet(I,'high')
printSet(I,'high',true)
printSet(I,'high',false)

% test normal set
a = [1;-1];
b = [2;3];
I = interval(a,b);

printSet(I)
printSet(I,'high')
printSet(I,'high',true)
printSet(I,'high',false)

% n-d arrays
lb = reshape([ 1.000 3.000 2.000 5.000 -3.000 0.000 2.000 1.000 0.000 -2.000 -1.000 3.000 0.000 0.000 0.000 0.000 1.000 -1.000 1.000 0.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
ub = reshape([ 1.500 4.000 4.000 10.000 -1.000 0.000 3.000 2.000 1.000 0.000 2.000 4.000 0.000 0.000 0.000 0.000 2.000 -0.500 3.000 2.000 0.000 0.000 0.000 0.000 ], [2,2,2,3]);
I = interval(lb,ub);
printSet(I)

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
