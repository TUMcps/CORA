function res = test_reachSet_isemptyobject
% test_reachSet_isemptyobject - unit test function for isemptyobject
%
% Syntax:
%    res = test_reachSet_isemptyobject()
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
% Written:       01-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty reachSet
R = reachSet();
res = isemptyobject(R);

% only time-point solution
timePoint.time = {0;0.01;0.02};
timePoint.set = {interval(0,1);interval(1,2);interval(2,3)};
R = reachSet(timePoint);
res(end+1,1) = ~isemptyobject(R);

% time-point and time-interval solution
timeInt.time = {interval(0,0.01),interval(0.01,0.02)};
timeInt.set = {interval(0,2),interval(1,3)};
R = reachSet(timePoint,timeInt);
res(end+1,1) = ~isemptyobject(R);

% 2 branches: one empty
R1 = reachSet();
R2 = reachSet(timePoint);
R = add(R1,R2);
res(end+1,1) = ~isemptyobject(R);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
