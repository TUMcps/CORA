function res = test_location_calcBasis
% test_location_calcBasis - test function for computation of orthogonal
%    basis for guard intersection
%
% Syntax:
%    res = test_location_calcBasis
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
% Written:       17-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% 2D location: flow moving to the right
inv = polytope([1,0],1);
guard = polytope([],[],[1,0],1);
reset = linearReset([1,0;0,1],[],[5;0]);
trans = transition(guard,reset,1);
flow = linearSys(zeros(2),0,[1;0]);
loc = location(inv,trans,flow);

% intersecting reachable sets
R{1} = zonotope([0;1],[1,0;0,1]);
R{2} = zonotope([1;0],[1,0;0,1]);
R{3} = zonotope([2;-1],[1,0;0,1]);

% go through different methods
options.enclose = {'box'};
B = calcBasis(loc,R,guard,options);
assert(all(all(withinTol(B{1},eye(2)))));

options.enclose = {'pca'};
B = calcBasis(loc,R,guard,options);
assert(all(all(withinTol(B{1},[0 1; 1 0]))));

options.enclose = {'flow'};
params.uTrans = [0;0];
params.w = 0;
B = calcBasis(loc,R,guard,options,params);
assert(all(all(withinTol(B{1},eye(2)))));

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
