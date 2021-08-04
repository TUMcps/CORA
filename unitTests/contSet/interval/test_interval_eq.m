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
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      29-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% TEST 1: Analytical ------------------------------------------------------
% create interval
lower1 = [-2; -4; -3];
upper1 = [ 2;  3;  1];
Int1 = interval(lower1, upper1);
lower2 = [-3;  0; -4];
upper2 = [ 2;  3;  1];
Int2 = interval(lower2, upper2);

% compute non-equality
res(1) = ~(Int1 == Int2); % different intervals
res(2) = Int1 == Int1;    % same interval
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
res_rand = Int1 == Int2;
res_rand_true = all(lower1 == lower2) && all(upper1 == upper2);
res(3) = res_rand == res_rand_true;
% -------------------------------------------------------------------------

% compare results
res = all(res);

if res
    disp('test_eq successful');
else
    disp('test_eq failed');
end

%------------- END OF CODE --------------