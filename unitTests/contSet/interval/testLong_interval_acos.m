function res = testLong_interval_acos
% testLong_interval_acos - unit test function of sine for intervals,
%    overloaded 'acos()' function for intervals
%
% Syntax:  
%    res = testLong_interval_acos
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
% See also: mtimes

% Author:       Mark Wetzlinger
% Written:      08-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% Random tests
numberRandTests = 10000;

% all x<1 or x>1 must return NaN
res_rand(1) = true;
for i=1:numberRandTests
    a = interval(-1-randi(1000,1,1)*rand(1),1+randi(1000,1,1)*rand(1));
    try
        b = acos(a); % error should be thrown here
        % in case no error is thrown, all should be NaN
        if ~all(isnan(b.inf)) || ~all(isnan(b.sup))
            res_rand(1) = false;
            break
        end
    catch
        continue;
    end
end

% all x>=-1 and x<=1 map to somewhere in [0,pi]
res_rand(2) = true;
for i=1:numberRandTests
    a = interval(-rand(1),rand(1));
    b = acos(a);
    if b.inf < 0 || b.sup > pi
        res_rand(2) = false;
        break
    end
end


% final test result
res = all(res_rand);

%------------- END OF CODE --------------