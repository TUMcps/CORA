function res = test_location_potOut
% test_location_potOut - test function for potOut
%
% Syntax:
%    res = test_location_potOut
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
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init location
inv = interval([-2;-1],[3;5]);
guard = conHyperplane([1 0],3);
reset = struct('A',eye(2),'c',[2;0]);
trans(1) = transition(guard,reset,2);
guard = conHyperplane([0 1],5);
reset = struct('A',eye(2),'c',[0;2]);
trans(2) = transition(guard,reset,2);
flow = linearSys(zeros(2),0,[1;1]);
loc = location(inv,trans,flow);

% cell-array of reachable sets
% -- not intersecting any guard
timePoint.set{1} = zonotope([1;0],0.05*eye(2));
timePoint.set{2,1} = zonotope([2;0],0.05*eye(2));
% -- only intersecting guard of first transition
timePoint.set{3} = zonotope([3;0],0.05*eye(2));
timePoint.set{4} = zonotope([3;2],0.05*eye(2));
timePoint.set{5} = zonotope([3;4],0.05*eye(2));
% -- intersecting both guards
timePoint.set{6} = zonotope([3;5],0.05*eye(2));
% -- only intersecting guard of second transition
timePoint.set{7} = zonotope([2;5],0.05*eye(2));
timePoint.set{8} = zonotope([1;5],0.05*eye(2));
% -- not intersecting any guard
timePoint.set{9} = zonotope([1;4],0.05*eye(2));
timePoint.time = {1;2;3;4;5;6;7;8;9};
for i=1:length(timePoint.set)-1
    timeInt.set{i,1} = convHull(timePoint.set{i},timePoint.set{i+1});
    timeInt.time{i,1} = {interval(timePoint.time{i},timePoint.time{i+1})};
end
R = reachSet(timePoint,timeInt);

% indices
minInd = [3;6];
maxInd = [6;8];

% evaluate function
R = potOut(loc,R,minInd,maxInd);

% sets 3 to 8 intersect and have been altered
Rtp = R.timePoint.set(3:8);
Rti = R.timeInterval.set(3:8);

% all intersecting sets must be polytopes
res = all(cellfun(@(x) isa(x,'polytope'),Rtp,'UniformOutput',true));
res(end+1,1) = all(cellfun(@(x) isa(x,'polytope'),Rti,'UniformOutput',true));
% all intersecting sets must be contained in invariant
res(end+1,1) = all(cellfun(@(x) contains(inv,x),Rtp,'UniformOutput',true));
res(end+1,1) = all(cellfun(@(x) contains(inv,x),Rti,'UniformOutput',true));


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
