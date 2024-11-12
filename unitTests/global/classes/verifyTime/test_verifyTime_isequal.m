function res = test_verifyTime_isequal
% test_verifyTime_isequal - unit test for helper class verifyTime
%
% Syntax:
%    res = test_verifyTime_isequal
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

% empty object
Itime_empty = verifyTime();
assert(isequal(Itime_empty,Itime_empty));

% single time interval
Itime = verifyTime([0,5]);
assert(isequal(Itime,Itime));
% ...with non-compact representation
Itime_alt = verifyTime([0,2;2,5]);
assert(isequal(Itime,Itime_alt));
assert(isequal(Itime_alt,Itime));
% ...different time interval
Itime_other = verifyTime([0,6]);
assert(~isequal(Itime,Itime_other));

% multiple disjoint time intervals
Itime1 = verifyTime([0,2; 2,5; 5,10; 12,15; 20,100]);
Itime2 = verifyTime([0,7; 7,10; 12,15; 20,100]);
assert(isequal(Itime1,Itime2));


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
