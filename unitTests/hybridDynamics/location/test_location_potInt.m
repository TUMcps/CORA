function res = test_location_potInt
% test_location_potInt - test function for finding potential intersections
%    of the reachable set and the guard set
%
% Syntax:
%    res = test_location_potInt
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

% reachSet object
% -- not intersecting any guard
sets.set{1} = zonotope([1;0],0.05*eye(2));
sets.set{2,1} = zonotope([2;0],0.05*eye(2));
% -- only intersecting guard of first transition
sets.set{3} = zonotope([3;0],0.05*eye(2));
sets.set{4} = zonotope([3;2],0.05*eye(2));
sets.set{5} = zonotope([3;4],0.05*eye(2));
% -- intersecting both guards
sets.set{6} = zonotope([3;5],0.05*eye(2));
% -- only intersecting guard of second transition
sets.set{7} = zonotope([2;5],0.05*eye(2));
sets.set{8} = zonotope([1;5],0.05*eye(2));
% -- not intersecting any guard
sets.set{9} = zonotope([1;4],0.05*eye(2));
% time is irrelevant
for i=1:9
    sets.time{i,1} = interval((i-1),i);
end
% init object (no time-point solutions)
R = reachSet([],sets);

% set options (dummy)
options.finalLoc = 4;

% check for potential intersections
[guards,setIndices] = potInt(loc,R,options);

% total intersections: 7
res = length(guards) == 7;
% correct guards for each intersection
res(end+1,1) = all(guards == [1,1,1,1,2,2,2]');
% correct indicies of intersecting sets
res(end+1,1) = all(setIndices == [3,4,5,6,6,7,8]');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
