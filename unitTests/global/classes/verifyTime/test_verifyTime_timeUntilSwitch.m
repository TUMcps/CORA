function res = test_verifyTime_timeUntilSwitch
% test_verifyTime_timeUntilSwitch - unit test for helper class verifyTime
%
% Syntax:
%    res = test_verifyTime_timeUntilSwitch
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false

% Authors:       Mark Wetzlinger
% Written:       10-November-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% tolerance
tol = 1e-12;

% empty object
Itime = verifyTime();
assert(isempty(timeUntilSwitch(Itime,2)));

% single time interval
Itime = verifyTime([2,5]);
tSwitch = timeUntilSwitch(Itime,1);
assert(withinTol(tSwitch,1,tol));
tSwitch = timeUntilSwitch(Itime,2);
assert(withinTol(tSwitch,3,tol));
tSwitch = timeUntilSwitch(Itime,3);
assert(withinTol(tSwitch,2,tol));
tSwitch = timeUntilSwitch(Itime,5);
assert(isempty(tSwitch));
tSwitch = timeUntilSwitch(Itime,6);
assert(isempty(tSwitch));

% multiple time intervals (non-minimal representation)
Itime = verifyTime([2,5; 10,20; 20,30]);
tSwitch = timeUntilSwitch(Itime,1);
assert(withinTol(tSwitch,1,tol));
tSwitch = timeUntilSwitch(Itime,2);
assert(withinTol(tSwitch,3,tol));
tSwitch = timeUntilSwitch(Itime,3);
assert(withinTol(tSwitch,2,tol));
tSwitch = timeUntilSwitch(Itime,5);
assert(withinTol(tSwitch,5,tol));
tSwitch = timeUntilSwitch(Itime,15);
assert(withinTol(tSwitch,15,tol));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
