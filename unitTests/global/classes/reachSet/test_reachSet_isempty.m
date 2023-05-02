function res = test_reachSet_isempty
% test_reachSet_isempty - unit test function for isempty
%
% Syntax:  
%    res = test_reachSet_isempty()
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

% Author:       Mark Wetzlinger
% Written:      01-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty reachSet
R = reachSet();
res = isempty(R);

% only time-point solution
timePoint.time = {0;0.01;0.02};
timePoint.set = {interval(0,1);interval(1,2);interval(2,3)};
R = reachSet(timePoint);
res(end+1,1) = ~isempty(R);

% time-point and time-interval solution
timeInt.time = {interval(0,0.01),interval(0.01,0.02)};
timeInt.set = {interval(0,2),interval(1,3)};
R = reachSet(timePoint,timeInt);
res(end+1,1) = ~isempty(R);

% 2 branches: one empty
R1 = reachSet();
R2 = reachSet(timePoint);
R = add(R1,R2);
res(end+1,1) = ~isempty(R);

% combine results
res = all(res);

%------------- END OF CODE --------------