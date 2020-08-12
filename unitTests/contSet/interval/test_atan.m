function res = test_atan
% test_atan - unit_test_function of sine for intervals - Overloaded 'atan()' function for intervals
%
% Syntax:  
%    res = test_atan
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


% 2. random tests ---------------------------------------------------------

numberRandTests = 10000;

% all values must be between -pi/2 and pi/2
res_rand(1) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    b = atan(a);
    if b.inf < -pi/2 || b.sup > pi/2
        res_rand(1) = false;
        break
    end
end

% all x<0 map to atan(x)<0, all x>0 map to atan(x)>0
res_rand(2) = true;
for i=1:numberRandTests
    a = interval(-randi(1000,1,1)*rand(1),randi(1000,1,1)*rand(1));
    b = atan(a);
    if b.inf >= 0 || b.sup <= 0
        res_rand(2) = false;
        break
    end
end



% final test result -------------------------------------------------------

res = all(res_val) && all(res_rand);

if res
    disp('test_atan successful');
else
    disp('test_atan failed');
end

end

%------------- END OF CODE --------------