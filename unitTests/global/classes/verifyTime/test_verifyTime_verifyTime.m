function res = test_verifyTime_verifyTime
% test_verifyTime_verifyTime - unit test for helper class verifyTime
%
% Syntax:
%    res = test_verifyTime_verifyTime
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

% no input arguments
Itime = verifyTime();
assert(all(size(Itime.bounds) == [0,2]));

% single time interval
bounds = [0,5];
Itime = verifyTime(bounds);
assert(all(Itime.bounds == bounds));

% multiple time intervals
bounds = [0,2; 3,4; 6,8; 10,20];
Itime = verifyTime(bounds);
assert(all(Itime.bounds == bounds,'all'));


% wrong initializations
% - wrong size
assertThrowsAs(@verifyTime,'CORA:wrongInputInConstructor',[0,1,2]);
% - non-finite values
assertThrowsAs(@verifyTime,'CORA:wrongInputInConstructor',[0,Inf]);
assertThrowsAs(@verifyTime,'CORA:wrongInputInConstructor',[0,1; NaN,2]);
% - lower bound is smaller than upper bound
assertThrowsAs(@verifyTime,'CORA:wrongInputInConstructor',[0,-1]);
assertThrowsAs(@verifyTime,'CORA:wrongInputInConstructor',[0,1; 3,2; 5,10]);
% - time intervals must be in order
assertThrowsAs(@verifyTime,'CORA:wrongInputInConstructor',[0,3; 2,3; 5,10]);


% test completed
res = true;

% ------------------------------ END OF CODE ------------------------------
