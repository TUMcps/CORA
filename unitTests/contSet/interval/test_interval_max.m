function res = test_interval_max
% test_interval_max - unit test function of dim
%
% Syntax:
%    res = test_interval_max
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Tobias Ladner
% Written:       16-December-2021
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% basic examples
I1 = interval([-2;-1],[2;1]);
Im = max(I1, 0);
res(end+1,1) = isequal(Im, interval([0;0], [2;1]));

I2 = interval([-1;1], [1;3]);
Im = max(I1, I2);
res(end+1,1) = isequal(Im, interval([-1;1], [2;3]));

% other set representations
Z = zonotope([1;-1], [2;1]);
Im = max(I1, Z);
res(end+1,1) = isequal(Im, interval([-1;-1], [3;1]));

% empty case
I = interval();
Im = max(I, I);
res(end+1,1) = representsa(Im,'emptySet');

I = interval();
Im = max(I, []);
res(end+1,1) = representsa(Im,'emptySet');

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
