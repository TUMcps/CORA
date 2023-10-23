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
guard = conHyperplane([1 0],4);
reset = struct('A',[2 -1; 1 3],'c',[1;-1]);
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
res = res_;
% reach set must only contain the first time-point solution
res(end+1,1) = ~isempty(R.timePoint) && length(R.timePoint.set) == 1;
% no time-interval solution
res(end+1,1) = isempty(R.timeInterval);
% only one jumped set
res(end+1,1) = length(Rjump) == 1;
% apply reset function to set
res(end+1,1) = isequal(Rjump{1}.set,reset.A*R0+reset.c);
% goal location = 2
res(end+1,1) = Rjump{1}.loc == trans.target;
% no time has passed
res(end+1,1) = Rjump{1}.time == tStart;

% combine results
res = all(res);


% ------------------------------ END OF CODE ------------------------------
