function res = test_reachSet_hasTimeInterval
% test_reachSet_hasTimeInterval - unit test function for hasTimeInterval
%
% Syntax:
%    res = test_reachSet_hasTimeInterval()
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

% Authors:       Tobias Ladner
% Written:       30-April-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% assume true
res = true;

% simple continuous-time linear system ------------------------------------
A = [0.1 1; -1 0.1];
B = 1;
sys = linearSys(A,B);

% model parameters and reachability options
params.tFinal = 1;
params.R0 = zonotope([10;5],0.2*eye(2));
options.linAlg = 'adaptive';
options.error = 0.1;

% compute reachable set
R = reach(sys,params,options);
assert(hasTimeInterval(R))

% simple discrete-time linear system --------------------------------------
A = [0.1 1; -1 0.1];
B = 1;
sys = linearSysDT(A,B,1);

% model parameters and reachability options
params.tFinal = 1;
params.R0 = zonotope([10;5],0.2*eye(2));
options = struct;
options.linAlg = 'adaptive';

% compute reachable set
R = reach(sys,params,options);
assert(~hasTimeInterval(R))

% special case: only initial set present ----------------------------------
% (happens if reachable set explodes in first step; ensures correct plotting)

timePoint = struct;
timePoint.set = {params.R0};
timePoint.time = {0};
timeInt = struct;
timeInt.set = cell(1,0);
timeInt.time = cell(1,0);
R = reachSet(timePoint,timeInt);
assert(~hasTimeInterval(R))

end


% ------------------------------ END OF CODE ------------------------------
