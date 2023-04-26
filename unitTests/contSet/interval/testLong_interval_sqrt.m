function res = testLong_interval_sqrt
% testLong_interval_sqrt - unit test function of sqrt
%
% Syntax:  
%    res = testLong_interval_sqrt
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

% Author:       Mark Wetzlinger
% Written:      29-August-2019
% Last update:  08-August-2020 (add random tests)
% Last revision:---

%------------- BEGIN CODE --------------

% random tests
numberRandTests = 5000;

% if inf < 0, result is NaN
res = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    try
        b = sqrt(a); % should throw an error
        % in case no error thrown, result should be NaN
        if ~all(b.inf <= 0)
            res = false; return;
        end
    catch
        continue
    end
end

% if inf > 0, result also inf > 0
for i=1:numberRandTests
    randinf = randi(1000,1,1)*rand(1);
    a = interval(randinf,randinf*2);
    b = sqrt(a);
    if any(b.inf <= 0)
        res = false;
        return;
    end
end

%------------- END OF CODE --------------