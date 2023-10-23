function res = test_reachSet_query
% test_reachSet_query - unit test function for query
%
% Syntax:
%    res = test_reachSet_query()
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
R = reachSet(timePoint,timeInt);

% queries: 'reachSet', 'reachSetTimePoint', 'finalSet', 'tVec'
R1 = query(R,'reachSet');
for i=1:length(R1)
    if ~isequal(R1{i},R.timeInterval.set{i})
        res = false; break
    end
end

R2 = query(R,'reachSetTimePoint');
for i=1:length(R2)
    if ~isequal(R2{i},R.timePoint.set{i})
        res = false; break
    end
end

Rfinal = query(R,'finalSet');
if ~isequal(Rfinal,R.timePoint.set{end})
    res = false;
end

tVec = query(R,'tVec');
for i=1:length(tVec)
    if ~withinTol(tVec(i),t(i+1)-t(i))
        res = false; break
    end
end

% ------------------------------ END OF CODE ------------------------------
