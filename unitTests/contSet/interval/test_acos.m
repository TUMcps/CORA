function res = test_acos
% test_acos - unit_test_function of sine for intervals - Overloaded 'acos()' function for intervals
%
% Syntax:  
%    res = test_acos
%
% Inputs:
%    no
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
% See also: mtimes

% Author:       Mark Wetzlinger
% Written:      08-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

% 1. Check special values -------------------------------------------------

% x     acos(x)
% -1    pi
% 0     pi/2
% 1     0

a = interval([-1,0],[1,0]);
b = acos(a);

binf_true = [0,pi/2];
bsup_true = [pi,pi/2];

res_val = true;
if any(abs( b.inf - binf_true ) > tol) || any(abs( b.sup - bsup_true ) > tol)
    res_val = false;
end


% 2. random tests ---------------------------------------------------------

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



% final test result -------------------------------------------------------

res = all(res_val) && all(res_rand);

if res
    disp('test_acos successful');
else
    disp('test_acos failed');
end

end

%------------- END OF CODE --------------