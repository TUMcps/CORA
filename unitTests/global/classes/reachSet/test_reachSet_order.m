function res = test_reachSet_order
% test_reachSet_order - unit test function for order
%
% Syntax:
%    res = test_reachSet_order()
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

% Authors:       Mark Wetzlinger
% Written:       10-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% instantiate set, propagation matrix
Z = zonotope([0;0],[1 0 -2; 2 -1 1]);
shift = [1;0];
steps = 10;
tPoint = [4;1;9;5;2;3;6;8;7;10];

% compute time-point solution
for i=1:steps
    timePoint.set{i,1} = Z + tPoint(i)*shift;
    timePoint.time{i,1} = tPoint(i);
end
R = reachSet(timePoint);

% mapped reach set
R_ordered = order(R);

% true solution
for i=1:steps
    timePoint_sort.set{i,1} = Z + i*shift;
    timePoint_sort.time{i,1} = i;
end
R_true = reachSet(timePoint_sort);

% compare results
res = isequal(R_true,R_ordered);

% ------------------------------ END OF CODE ------------------------------
