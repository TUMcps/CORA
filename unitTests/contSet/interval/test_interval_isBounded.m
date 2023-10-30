function res = test_interval_isBounded
% test_interval_isBounded - unit test function of isBounded
%
% Syntax:
%    res = test_interval_isBounded
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

% Authors:       Tobias Ladner
% Written:       26-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% instantiate bounded intervals
I = interval();
resvec(end+1) = isBounded(I);

I = interval(1);
resvec(end+1) = isBounded(I);

I = interval([-2;-1;-3],[1;1;2]);
resvec(end+1) = isBounded(I);

% instantiate unbounded intervals
I = interval(inf);
resvec(end+1) = ~isBounded(I);

I = interval([-2;-inf;-3],[1;1;2]);
resvec(end+1) = ~isBounded(I);

I = interval([-2;-inf;-3],[1;inf;2]);
resvec(end+1) = ~isBounded(I);

% check results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
