function res = test_reachSet_find
% test_reachSet_find - unit test function for find
%
% Syntax:
%    res = test_reachSet_find()
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

% init result
res = true;

% instantiate set, propagation matrix
Z = zonotope([1;1;4],[1 0 -2; 2 -1 1; -2 4 -1]);
s = [1; -1; 2];
t = [0;2;5;10;13;17;24;33;39;40;50];
parents = [0, 1];
locs = [1, 2];
steps = 10;

timePoint.set{1,1} = Z;
timePoint.time{1,1} = t(1);
% propagate sets
for i=1:steps
    timePoint.set{i+1,1} = Z + i*s;
    timePoint.time{i+1,1} = t(i+1);
    timeInt.set{i,1} = Z - i*s;
    timeInt.time{i,1} = interval(t(i),t(i+1));
end
R1 = reachSet(timePoint,timeInt,parents(1),locs(1));
R2 = R1 + s;
R2 = reachSet(R2.timePoint,R2.timeInterval,parents(2),locs(2));
R = [R1;R2];

% properties: 'location', 'parent', 'time'
R_loc = find(R,'location',locs(2));
if ~isequal(R_loc,R2)
    res = false;
end

R_parent = find(R,'parent',parents(1));
if ~isequal(R_parent,R1)
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
