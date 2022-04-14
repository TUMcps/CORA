function res = test_interval_acos
% test_interval_acos - unit test function of sine for intervals,
%    overloaded 'acos()' function for intervals
%
% Syntax:  
%    res = test_interval_acos
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


% final test result -------------------------------------------------------

res = all(res_val);

if res
    disp('test_acos successful');
else
    disp('test_acos failed');
end

end

%------------- END OF CODE --------------