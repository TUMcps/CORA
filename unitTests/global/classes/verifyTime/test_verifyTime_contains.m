function res = test_verifyTime_contains
% test_verifyTime_contains - unit test for helper class verifyTime
%
% Syntax:
%    res = test_verifyTime_contains
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
Itime = verifyTime();
t = 2;
assert(~contains(Itime,t));

% single time interval
Itime = verifyTime([2,5]);
[isContained,idx] = contains(Itime,1);
assert(~isContained && isempty(idx));
[isContained,idx] = contains(Itime,3);
assert(isContained && idx == 1);
[isContained,idx] = contains(Itime,7);
assert(~isContained && isempty(idx));

% multiple time intervals (non-minimal representation)
Itime = verifyTime([2,5; 10,20; 20,30]);
[isContained,idx] = contains(Itime,1);
assert(~isContained && isempty(idx));
[isContained,idx] = contains(Itime,3);
assert(isContained && idx == 1);
[isContained,idx] = contains(Itime,20);
% note: we return first interval, for which containment holds true
assert(isContained && idx == 2);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
