function res = test_verifyTime_compact
% test_verifyTime_compact - unit test for helper class verifyTime
%
% Syntax:
%    res = test_verifyTime_compact
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
Itime_compact = compact(Itime);
assert(isequal(Itime,Itime_compact));

% single time interval
bounds = [0,5];
Itime = verifyTime(bounds);
Itime_compact = compact(Itime);
bounds_true = bounds;
assert(compareMatrices(Itime_compact.bounds,bounds_true,tol,"equal",true));

% multiple time intervals -> single time interval
bounds = [0,5; 5,10; 10,20];
Itime = verifyTime(bounds);
Itime_compact = compact(Itime);
bounds_true = [0,20];
assert(compareMatrices(Itime_compact.bounds,bounds_true,tol,"equal",true));

% multiple time intervals -> no simplification
bounds = [0,5; 6,10; 12,20];
Itime = verifyTime(bounds);
Itime_compact = compact(Itime);
bounds_true = bounds;
assert(compareMatrices(Itime_compact.bounds,bounds_true,tol,"equal",true));

% multiple time intervals -> partial simplification
bounds = [0,2; 2,5; 5,10; 12,15; 20,30; 30,50];
Itime = verifyTime(bounds);
Itime_compact = compact(Itime);
bounds_true = [0,10; 12,15; 20,50];
assert(compareMatrices(Itime_compact.bounds,bounds_true,tol,"equal",true));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
