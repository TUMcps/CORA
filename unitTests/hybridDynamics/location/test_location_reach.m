function res = test_location_reach
% test_location_reach - test function for reach
%
% Syntax:
%    res = test_location_reach
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
inv = polytope([1 0],2);
guard = conHyperplane([1 0],2);
reset = struct('A',eye(2),'c',[3;-2]);
trans = transition(guard,reset,2);
flow = linearSys(zeros(2),0,[1;-0.2]);
loc = location(inv,trans,flow);

% start set and time
R0 = zonotope([0.5;3],0.5*eye(2));
tStart = interval(2,2.5);

% set options
options.R0 = R0;
options.timeStep = 1;
options.taylorTerms = 2;
options.zonotopeOrder = 10;
options.U = zonotope(zeros(2,1));
options.tStart = 2;
options.tFinal = 10;
options.specification = specification();
options.finalLoc = 4;
options.guardIntersect = 'polytope';
options.enclose = {'box'};
options.intersectInvariant = false;

% call reach
[R,Rjump,res_] = reach(loc,R0,tStart,options);

% no specification violated
res = res_;
% reachable set needs to consist of four sets
res(end+1,1) = length(R.timeInterval.set) == 4;
% all Rjump - reset.c need to intersect the guard set
res(end+1,1) = arrayfun(@(x) isIntersecting(x.set-reset.c,guard),Rjump,...
    'UniformOutput',true);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
