function res = test_location_instantReset
% test_location_instantReset - test function for instantReset
%
% Syntax:
%    res = test_location_instantReset
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
inv = interval([-2;1],[1;4]);
guard = polytope([],[],[1 0],4);
reset = linearReset([2 -1; 1 3],[],[1;-1]);
trans = transition(guard,reset,2);
flow = linearSys([2 1; -1 2],1);
loc = location(inv,trans,flow);

% start set and time
R0 = zonotope([4;2]); 
tStart = 2;

% options
options.instantTransition = true;
options.U = zonotope(0);

% reset
[R,Rjump,res_] = instantReset(loc,R0,tStart,options);

% 'res' must be true as no specification was violated
assert(res_);
% reach set must only contain the first time-point solution
assert(~isempty(R.timePoint) && length(R.timePoint.set) == 1);
% no time-interval solution
assert(isempty(R.timeInterval));
% only one jumped set
assert(length(Rjump) == 1);
% apply reset function to set
assert(isequal(Rjump{1}.set,reset.A*R0+reset.c));
% goal location = 2
assert(Rjump{1}.loc == trans.target);
% no time has passed
assert(Rjump{1}.time == tStart);

% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
