function res = test_interval_enlarge
% test_interval_enlarge - unit test function of enlarge
%
% Syntax:  
%    res = test_interval_enlarge
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
lower = [-2; -4; -3];
upper = [ 2;  3;  1];
factor = 2;
Int = interval(lower, upper);

% compute enlarge
IntEnlarge = enlarge(Int, factor);

% true size
lower_true = [-4; -7.5; -5];
upper_true = [ 4;  6.5;  3];
IntEnlarge_true = interval(lower_true, upper_true);

% compare results
res_analytical = all(infimum(IntEnlarge) == infimum(IntEnlarge_true)) && ...
        all(supremum(IntEnlarge) == supremum(IntEnlarge_true));
% -------------------------------------------------------------------------

% TEST 2: Random ----------------------------------------------------------
% create random interval
dim = floor(1 + 9*rand(1));
lower = -10*rand(dim,1);
upper = 10*rand(dim,1);
factor = 2;
Int = interval(lower, upper);

% compute enlarge
IntEnlarge = enlarge(Int, factor);

% true size
halfWay = (upper + lower)/2;
edgeLength = (upper - lower)/2;
lower_true = halfWay - edgeLength*factor;
upper_true = halfWay + edgeLength*factor;
IntEnlarge_true = interval(lower_true, upper_true);

% compare results
res_rand = all(infimum(IntEnlarge) == infimum(IntEnlarge_true)) && ...
        all(supremum(IntEnlarge) == supremum(IntEnlarge_true));
% -------------------------------------------------------------------------

% add results
res = res_analytical && res_rand;

if res
    disp('test_enlarge successful');
else
    disp('test_enlarge failed');
end

%------------- END OF CODE --------------