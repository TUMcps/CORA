function res = testLongDuration_interval_sqrt
% testLongDuration_interval_sqrt - unit test function of sqrt
%
% Syntax:  
%    res = testLongDuration_interval_sqrt
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

% Random tests ---------------------------------------------------------

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
res = all(res_rand);

if res
    disp('test_sqrt successful');
else
    disp('test_sqrt failed');
end

%------------- END OF CODE --------------