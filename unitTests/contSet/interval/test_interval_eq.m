function res = test_interval_eq
% test_interval_eq - unit test function of eq, overloaded '==' operator
%
% Syntax:
%    res = test_interval_eq
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
% See also: -

% Authors:       Mark Wetzlinger
% Written:       29-August-2019
% Last update:   04-December-2023 (MW, add empty and unbounded cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
I1 = interval.empty(2);
res(end+1,1) = I1 == I1;


% bounded
I1 = interval([-2; -4; -3],[2; 3; 1]);
I2 = interval([-3; 0; -4],[2; 3; 1]);
res(end+1,1) = I1 == I1;
res(end+1,1) = ~(I1 == I2);

% unbounded
I1 = interval(-Inf,0);
I2 = interval(-Inf,Inf);
I3 = interval(0,Inf);
res(end+1,1) = I1 == I1;
res(end+1,1) = I2 == I2;
res(end+1,1) = I3 == I3;
res(end+1,1) = ~(I1 == I2);
res(end+1,1) = ~(I1 == I3);
res(end+1,1) = ~(I2 == I3);


% compare results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
