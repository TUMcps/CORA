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
guard = polytope([],[],[1 0],2);
reset = linearReset(eye(2),zeros(2),[3;-2]);
trans = transition(guard,reset,2);
flow = linearSys(zeros(2),zeros(2),[1;-0.2]);
loc = location(inv,trans,flow);

% start set and time
R0 = zonotope([0.5;3],0.5*eye(2));
tStart = interval(2,2.5);

% set parameters and options
params.R0 = R0;
params.U = zonotope(zeros(2,1));
params.tStart = tStart;
params.tFinal = 10;
params.finalLoc = 4;

options.timeStep = 1;
options.taylorTerms = 2;
options.zonotopeOrder = 10;
options.specification = specification();
options.guardIntersect = 'polytope';
options.enclose = {'box'};
options.intersectInvariant = false;

% call reach
[R,Rjump,res_] = reach(loc,params,options);

% no specification violated
assert(res_);
% reachable set needs to consist of four sets
assert(length(R.timeInterval.set) == 4);
% all Rjump - reset.c need to intersect the guard set
assert(arrayfun(@(x) isIntersecting(x.set-reset.c,guard),Rjump,...
    'UniformOutput',true));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
