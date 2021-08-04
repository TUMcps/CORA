function res = test_interval_atan
% test_interval_atan - unit test function of sine for intervals,
%    overloaded 'atan()' function for intervals
%
% Syntax:  
%    res = test_interval_atan
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

% x                         atan(x)     ..in deg
% 0                         0           0
% sqrt(3)/3                 pi/6        30
% 1                         pi/4        45
% sqrt(3)                   pi/3        60
% 2 - sqrt(3)               pi/12       15
% 2 + sqrt(3)               5pi/12      75
% sqrt(5-2*sqrt(5))         pi/5        36
% sqrt(5+2*sqrt(5))         2pi/5       72

a = interval([0,sqrt(3)/3],[1,sqrt(3)]);
b = atan(a);

binf_true = [0,pi/6];
bsup_true = [pi/4,pi/3];

res_val(1) = true;
if any(abs( b.inf - binf_true ) > tol) || any(abs( b.sup - bsup_true ) > tol)
    res_val(1) = false;
end

a = interval([2-sqrt(3);sqrt(5-2*sqrt(5))],[2+sqrt(3);sqrt(5+2*sqrt(5))]);
b = atan(a);

binf_true = [pi/12;pi/5];
bsup_true = [5*pi/12;2*pi/5];

res_val(2) = true;
if any(abs( b.inf - binf_true ) > tol) || any(abs( b.sup - bsup_true ) > tol)
    res_val(2) = false;
end



% final test result -------------------------------------------------------

res = all(res_val);

if res
    disp('test_atan successful');
else
    disp('test_atan failed');
end

end

%------------- END OF CODE --------------