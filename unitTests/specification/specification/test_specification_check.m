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
set = halfspace([1 1]/sqrt(2),1);
specUnsafe = specification(set,'unsafeSet');
specSafe = specification(set,'safeSet');
specInvariant = specification(set,'invariant');


% check set
Z = zonotope([-10;-8],[1 0 -2; 2 -1 1]);
res = ~check(specUnsafe,Z);
res(end+1,1) = check(specSafe,Z);
res(end+1,1) = check(specInvariant,Z);

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
res(end+1,1) = ~check(specUnsafe,R);
res(end+1,1) = check(specSafe,R);
res(end+1,1) = check(specInvariant,R);

% check simulation results
res(end+1,1) = ~check(specUnsafe,simRes);
res(end+1,1) = check(specSafe,simRes);
res(end+1,1) = check(specInvariant,simRes);


% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
