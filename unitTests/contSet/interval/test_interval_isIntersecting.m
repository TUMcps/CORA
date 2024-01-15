function res = test_interval_isIntersecting
% test_interval_isIntersecting - unit test function of isIntersecting
%    note: only interval-to-interval tested!
%
% Syntax:
%    res = test_interval_isIntersecting
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
% Written:       12-March-2021
% Last update:   04-December-2023 (MW, add more cases)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% empty case
I1 = interval.empty(1);
I2 = interval(-1,1);
res(end+1,1) = ~isIntersecting(I1,I2);
res(end+1,1) = ~isIntersecting(I2,I1);

% bounded
I1 = interval([-2;-1],[1;2]);
I2 = interval([-4;-2],[-1;0]);
res(end+1,1) = isIntersecting(I1,I2);
res(end+1,1) = isIntersecting(I2,I1);
I1 = interval([-2;-1],[1;2]);
I2 = interval([-4;-2],[-3;0]);
res(end+1,1) = ~isIntersecting(I1,I2);
res(end+1,1) = ~isIntersecting(I2,I1);

% bounded and unbounded
I1 = interval([-2;-1],[1;2]);
I2 = interval([-4;1],[-1;Inf]);
res(end+1,1) = isIntersecting(I1,I2);
res(end+1,1) = isIntersecting(I2,I1);

% unbounded and unbounded
I1 = interval(-Inf,0);
I2 = interval(-1,Inf);
res(end+1,1) = isIntersecting(I1,I2);
res(end+1,1) = isIntersecting(I2,I1);


% combine results
res = all(res);


% dimension mismatch
I1 = interval(-1,1);
I2 = interval([-1;-2],[2;1]);
try
    isIntersecting(I1,I2);
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
