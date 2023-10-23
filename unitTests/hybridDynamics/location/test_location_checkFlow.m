function res = test_location_checkFlow
% test_location_checkFlow - test function for checking the flow of
%    intersecting reachable sets to determine which ones can be exempt from
%    reset function (to avoid infinite splitting due to over-approximation)
%
% Syntax:
%    res = test_location_checkFlow
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

% init location with simple dynamics
inv = polytope([1 0],0);
guard = conHyperplane([1 0],0);
reset = struct('A',eye(2),'c',zeros(2,1));
trans = transition(guard,reset,1);
flow = linearSys(zeros(2),0,[-1;0]);
loc = location(inv,trans,flow);

% intersecting sets
c = [-1;1];
G = [0.4 -0.3 0 -0.5; 0.4 0 -0.2 0.6];
R{1} = zonotope(c,G);
R{2} = R{1} + [1;-0.7];
R{3} = R{2} + [1;-0.7];

% set required options
options.U = zonotope(0);

% check flow
[res_,R_] = checkFlow(loc,guard,R,options);
% no sets should be left
res = ~res_ && isempty(R_);

% init location with level set
syms x y;
eq = y - x^2;
inv = levelSet(eq,[x;y],'<=');
guard = levelSet(eq,[x;y],'==');
reset = struct('A',eye(2),'c',zeros(2,1));
trans = transition(guard,reset,1);
flow = linearSys(zeros(2),0,[1;1]);
loc = location(inv,trans,flow);

% intersecting sets
R{1} = interval([-2;-1],[-1;1]);
R{2} = interval([-0.5;-1],[0.5;1]);
R{3} = interval([1;-1],[2;1]);

% check flow
[res_,R_] = checkFlow(loc,guard,R,options);
% all but last set should be left
res(end+1,1) = res_ && length(R_) == 2 && ...
    isequal(R_{1},R{1}) && isequal(R_{2},R{2});

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
