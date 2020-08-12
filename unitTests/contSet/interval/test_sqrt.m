function res = test_sqrt
% test_sqrt - unit test function of sqrt
%
% Syntax:  
%    res = test_sqrt
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

% 2. random tests ---------------------------------------------------------

numberRandTests = 5000;

% if inf < 0, result is NaN
res_rand(1) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    try
        b = sqrt(a); % should throw an error
        % in case no error thrown, result should be NaN
        if ~all(b.inf <= 0)
            res_rand(1) = false;
            return;
        end
    catch
        continue
    end
end

% if inf > 0, result also inf > 0
res_rand(1) = true;
for i=1:numberRandTests
    randinf = randi(1000,1,1)*rand(1);
    a = interval(randinf,randinf*2);
    b = sqrt(a);
    if any(b.inf <= 0)
        res_rand(2) = false;
        return;
    end
end


% add results
res = all(res_val) && all(res_rand);

if res
    disp('test_sqrt successful');
else
    disp('test_sqrt failed');
end

%------------- END OF CODE --------------