function res = test_linearSys_verify_01_temporalLogic
% test_linearSys_verify_01_temporalLogic - unit test for automated 
%    verification of signal temporal logic specifications
%
% Syntax:
%    res = test_linearSys_verify_01_temporalLogic
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Niklas Kochdumper
% Written:       22-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = false;

% Analytical Test (Until) -------------------------------------------------

% system dynamics
A = [0 -1; 1 0];
B = [0;0];

sys = linearSys(A,B);

% reachability parameter
params.R0 = zonotope([0;-1],diag([0.1,0.1]));
params.U = zonotope(0);
params.tFinal = 2;

options = struct;
options.verifyAlg = 'stl:kochdumper';

% safe specification
x = stl('x',2);
eq = until(x(2) < -0.5,x(1) > 0.5,interval(0,2));

specSafe = specification(eq,'logic');

% unsafe specification
x = stl('x',2);
eq = until(x(2) < -0.7,x(1) > 0.7,interval(0,2));

specUnsafe = specification(eq,'logic');

% automated verification
resSafe = verify(sys,params,options,specSafe);
resUnsafe = verify(sys,params,options,specUnsafe);

% check results
if ~resSafe || resUnsafe
    error('Analytical test for the "until" operator failed!');
end


% Analytical Test (Finally) -----------------------------------------------

% system dynamics
A = [0 -1; 1 0];
B = [0;0];

sys = linearSys(A,B);

% reachability parameter
params.R0 = zonotope([0;-1],diag([0.1,0.1]));
params.U = zonotope(0);
params.tFinal = 2;

options = struct;
options.verifyAlg = 'stl:kochdumper';

% safe specification
x = stl('x',2);
goal = 0.2*interval(-[1;0.1],[1;0.1]) + [1;0];
eq = finally(in(x,goal),interval(0,2));

specSafe = specification(eq,'logic');

% unsafe specification
x = stl('x',2);
goal = 0.2*interval(-[1;0.1],[1;0.1]) + [1.2;0];
eq = finally(in(x,goal),interval(0,2));

specUnsafe = specification(eq,'logic');

% automated verification
resSafe = verify(sys,params,options,specSafe);
resUnsafe = verify(sys,params,options,specUnsafe);

% check results
if ~resSafe || resUnsafe
    error('Analytical test for the "finally" operator failed!');
end

% all checks ok
res = true;

end

% ------------------------------ END OF CODE ------------------------------
