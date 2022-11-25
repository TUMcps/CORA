function res = test_interval_sqrt
% test_interval_sqrt - unit test function of sqrt
%
% Syntax:  
%    res = test_interval_sqrt
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
% Last update:  08-August-2020 (add random tests)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% 1. analytical tests -----------------------------------------------------

% create interval (all positive)
lower = [4; 9; 4; 16; 1];
upper = [9; 25; 36; 100; 4];
Int = interval(lower, upper);

% compute sqrt
Int_sqrt = sqrt(Int);

% true solution
lower_true = [2; 3; 2; 4; 1];
upper_true = [3; 5; 6; 10; 2];
Int_true = interval(lower_true, upper_true);

% compare results
res_val(1) = Int_sqrt == Int_true;


% values of inf / sup below zero should throw error
dim = 5;
lower_rand = -rand(dim,1);
upper_rand = rand(dim,1);
Int = interval(lower_rand, upper_rand);
try
    Int_sqrt = sqrt(Int);
    res_val(2) = false;
catch
    res_val(2) = true;
end


% add results
res = all(res_val);

if res
    disp('test_sqrt successful');
else
    disp('test_sqrt failed');
end

%------------- END OF CODE --------------