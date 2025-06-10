function res = test_specification_check
% test_specification_check - unit test for check
%
% Syntax:
%    res = test_specification_check
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
% Written:       30-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init specifications
set = polytope([1 1]/sqrt(2),1);
specUnsafe = specification(set,'unsafeSet');
specSafe = specification(set,'safeSet');
specInvariant = specification(set,'invariant');


% check set
Z = zonotope([-10;-8],[1 0 -2; 2 -1 1]);
assert(~check(specUnsafe,Z));
assert(check(specSafe,Z));
assert(check(specInvariant,Z));

% init system
sys = linearSys([-4 -1; 1 -4],1);
params.R0 = zonotope([-100;-100],eye(2));
params.tFinal = 2;

% reachability 
options.timeStep = 0.02;
options.taylorTerms = 5;
options.zonotopeOrder = 20;
R = reach(sys,params,options);

% simulation
simOpt.points = 10;
simRes = simulateRandom(sys,params,simOpt);

% check reachable set
assert(~check(specUnsafe,R));
assert(check(specSafe,R));
assert(check(specInvariant,R));

% check simulation results
assert(~check(specUnsafe,simRes));
assert(check(specSafe,simRes));
assert(check(specInvariant,simRes));

% check timed ---

specTimed = specification(interval(-1,1),'safeSet',interval(1));

% no simulation in time, something is wrong with simRes object
simRes = simResult({ones(11,1)},{(0:0.01:0.1)'});
assert(~check(specTimed,simRes))
% but with time, all is good
simRes = simResult({ones(11,1)},{(0:0.1:1)'});
assert(check(specTimed,simRes))

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
