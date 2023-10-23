function res = test_interval_ne
% test_interval_ne - unit test function of ne, overloaded '~=' operator
%
% Syntax:
%    res = test_interval_ne
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
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% TEST 1: Analytical ------------------------------------------------------
% create interval
lower1 = [-2; -4; -3];
upper1 = [ 2;  3;  1];
Int1 = interval(lower1, upper1);
lower2 = [-3;  0; -4];
upper2 = [ 2;  3;  1];
Int2 = interval(lower2, upper2);

% compute non-equality
res(1) = Int1 ~= Int2;    % different intervals
res(2) = ~(Int1 ~= Int1); % same interval
% -------------------------------------------------------------------------

% TEST 2: Random ----------------------------------------------------------
% create interval
dim = floor(1 + 9*rand(1));
lower1 = -10*rand(dim,1);
upper1 = 10*rand(dim,1);
Int1 = interval(lower1, upper1);
lower2 = -10*rand(dim,1);
upper2 = 10*rand(dim,1);
Int2 = interval(lower2, upper2);

% compute non-equality
res_rand = Int1 ~= Int2;
res_rand_true = any(lower1 ~= lower2) || any(upper1 ~= upper2);
res(3) = res_rand == res_rand_true;
% -------------------------------------------------------------------------

% compare results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
