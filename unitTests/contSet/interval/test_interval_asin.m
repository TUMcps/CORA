function res = test_interval_asin
% test_interval_asin - unit test function of sine for intervals,
%    overloaded 'asin()' function for intervals
%
% Syntax:  
%    res = test_interval_asin
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
% See also: mtimes

% Author:       Mark Wetzlinger
% Written:      08-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

tol = 1e-9;
res = true;

% 1. Check special values -------------------------------------------------

% x     asin(x)
% -1    -pi/2
% 0     0
% 1     pi/2

a = interval([-1,0],[1,0]);
b = asin(a);

binf_true = [-pi/2,0];
bsup_true = [pi/2,0];

res_val = true;
if any(abs( b.inf - binf_true ) > tol) || any(abs( b.sup - bsup_true ) > tol)
    res_val = false;
end



% final test result -------------------------------------------------------

res = all(res_val);

if res
    disp('test_asin successful');
else
    disp('test_asin failed');
end

end

%------------- END OF CODE --------------