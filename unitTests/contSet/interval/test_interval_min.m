function res = test_interval_min
% test_interval_min - unit test function of dim
%
% Syntax:
%    res = test_interval_min
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

res = true;

% basic examples
I1 = interval([-2;-1],[2;1]);
Im = min(I1, 0);
res = res && isequal(Im, interval([-2;-1], [0;0]));

I2 = interval([-1;1], [1;3]);
Im = min(I1, I2);
res = res && isequal(Im, interval([-2;-1], [1;1]));

% other set representations
Z = zonotope([1;-1], [2;1]);
Im = min(I1, Z);
res = res && isequal(Im, interval([-2;-2], [2;0]));

% empty case
I = interval();
Im = min(I, I);
res = res & representsa_(Im,'emptySet',eps);

I = interval();
Im = min(I, []);
res = res & representsa_(Im,'emptySet',eps);

% ------------------------------ END OF CODE ------------------------------
